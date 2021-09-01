__doc__ = """
Takes in a Desmond MD trajectory, splits it into individual frame snapshots,
and runs each one through MMGBSA (after deleting waters and separating the
ligand from the receptor).

Copyright Schrodinger, LLC. All rights reserved.

edit: add POPC and other lipids to molecules that need to be deleted
"""

import os
import sys
import numpy
import argparse
import csv
from collections import OrderedDict, defaultdict
from schrodinger import structure
from schrodinger.utils import subprocess
from schrodinger.structutils import analyze
from schrodinger.job import jobcontrol, queue
from schrodinger.trajectory.desmondsimulation import create_simulation

RUN = os.path.join(os.environ['SCHRODINGER'], 'run')
PRIME_MMGBSA = os.path.join(os.environ['SCHRODINGER'], 'prime_mmgbsa')
MUT_PROP = "s_bioluminate_Mutations"
MMGBSA_PROPS = OrderedDict([
    ("dG", "r_psp_MMGBSA_dG_Bind"),
    ("dG(NS)", "r_psp_MMGBSA_dG_Bind(NS)"),
])
RESIDUE_SCANNING_PROPS = OrderedDict([
    ("dAffinity", "r_bioluminate_delta_Affinity"),
    ("dStability", "r_bioluminate_delta_Stability"),
    ("dSASA(total)", "r_bioluminate_delta_SASA_(total)"),
    ("dSASA(nonpolar)", "r_bioluminate_delta_SASA_(nonpolar)"),
    ("dSASA(polar)", "r_bioluminate_delta_SASA_(polar)"),
])
# FIXME include more MMGBSA and ResidueScanning properties?


def find_ligand(struct):
    """
    Detect the ligand in the given complex CT, and return an ASL for it.
    """

    def compare_natom(lig1, lig2):
        return cmp(len(lig1.atom_indexes), len(lig2.atom_indexes))

    ligand_searcher = analyze.AslLigandSearcher()
    ligand_searcher.min_atom_count = 15
    ligand_searcher.max_atom_count = 300
    ligand_searcher.exclude_ions = True
    ligand_searcher.exclude_amino_acids = False
    ligand_list = ligand_searcher.search(struct)

    if not ligand_list:
        print "ERROR: No ligand was detected in the input trajectory"
        sys.exit(1)

    if len(ligand_list) == 1:
        ligand = ligand_list[0]
        print "Found ligand (%i atoms)" % len(ligand.atom_indexes)
    else:
        # Select the ligand with the most atoms:
        print "Found %i potential ligands" % len(ligand_list)
        ligand_list.sort(compare_natom)
        ligand = ligand_list[0]
        print "Select the one with most atoms (%i)" % len(ligand.atom_indexes)

    return ligand.ligand_asl


def calc_ave_and_std(values):
    """
    For a given list of values, calculate and return the average and the
    standard deviation.
    """

    array = numpy.array(values)
    ave = numpy.mean(array, axis=0)
    std = numpy.std(array, axis=0)
    return ave, std


def run_residue_scanning(infile, res_file, lig_asl, basename):
    """
    Run Residue Scanning jobs, in parallel, on the input complex structures.

    @type infile: str
    @param infile: Path to the input structure file (multiple complexes).

    @type res_file: str
    @param res_file: Path to the complex input file.

    @type lig_asl: str
    @param lig_asl: ASL defining the ligand, for affinity calculation.

    @type basename: str
    @param basename: Base job name to use.

    @rtype: str
    @return: Path to the output file.
    """

    # Will use the host specified via the -HOST argument:
    jdj = queue.JobDJ(verbosity="verbose")

    out_files = []
    for i, st in enumerate(structure.StructureReader(infile), start=1):
        subjobname = '%s-rs-%03d' % (basename, i)
        in_file = '%s-in.maegz' % subjobname
        out_file = '%s-out.maegz' % subjobname
        st.write(in_file)
        out_files.append(out_file)
        cmd = [RUN, "residue_scanning_backend.py", "-jobname", subjobname,
            "-calc", "sasa_polar,sasa_nonpolar,sasa_total",
            "-fast", "-refine_mut", "prime_residue",
            "-receptor_asl", "NOT (%s)" % lig_asl,
            "-res_file", res_file, in_file]
        jdj.addJob(cmd)

    # Run all Residue Scanning jobs in parallel:
    jdj.run()

    merged_out_file = basename + "-residue-scanning-out.maegz"
    writer = structure.StructureWriter(merged_out_file)
    for out_file in out_files:
        for i, st in enumerate(structure.StructureReader(out_file)):
            if i == 0:
                continue # Skip the original structures
            writer.append(st)
    writer.close()
    return merged_out_file

if __name__ == '__main__':

    script_desc = """
        $SCHRODINGER/run thermal_mmgbsa.py <trajectory.cms>

        Specify the -HOST option to run on another host or to use multiple CPUs.
    """
    parser = argparse.ArgumentParser(description=script_desc)
    parser.add_argument("trajectory",
                        help="Input trajectory")

    parser.add_argument("-lig_asl",
                        help="ASL for the ligand. If not specified, ligand will be detected automatically.",
                        default=None)

    parser.add_argument("-j", "-o",
                        dest='jobname',
                        default='trajectory-mmgbsa',
                        help="Job name and base name for the output files. "
                        "Default is 'trajectory-mmgbsa'.")
    parser.add_argument("-start_frame", default=0, type=int,
                        help="Frame to start analysis from")
    parser.add_argument("-end_frame", type=int,
                        help="Frame to end analysis on")
    parser.add_argument("-step_size", default=1, type=int,
                        help="Analyze every Nth frame (default is 1)")
    parser.add_argument("-res_file",
                        help="Perform these mutations (Residue Scanning) before "
                        "running prime_mmgbsa. For format details, see: "
                        "$SCHRODINGER/run residue_scanning_backend.py -h")

    # The -HOST option actually gets parsed by $SCHRODINGER/run
    parser.add_argument("-HOST",
                        dest='ignored',
                        default='localhost',
                        help="Host to run the Prime jobs on. (e.g. 'localhost:4')")

    parser.add_argument("-NJOBS",
                        dest='NJOBS',
                        type=int,
                        help="Number of Prime subjobs to create.")

    try:
        args = parser.parse_args()
    except IOError as err:
        parser.error(str(err))
        sys.exit(1)

    if not args.trajectory:
        parser.error("Please specify input trajectory")
        sys.exit(1)

    in_cms = args.trajectory

    if in_cms.endswith('.tgz'):
        print "ERROR: Input trajectory must be uncompressed"
        sys.exit(1)
    # FIXME add ability to decompress trajectories

    if not in_cms.endswith('-out.cms'):
        print "ERROR: Unrecognized trajectory format. Must end with '.out.cms'."
        sys.exit(1)

    in_basename = in_cms[:-8]
    in_trj = in_basename + '_trj'

    if not os.path.isdir(in_trj):
        print "ERROR: Missing directory: '%s'" % in_trj
        sys.exit(1)

    csim = create_simulation(in_cms, in_trj)

    complexes_file = args.jobname + "-complexes.maegz"
    writer = structure.StructureWriter(complexes_file)

    water_atoms = None

    start_frame = args.start_frame
    end_frame = args.end_frame or csim.total_frame
    print 'Will process frames %i-%i (step size: %i)' \
          % (start_frame, end_frame, args.step_size)

    if args.lig_asl:
        print "Will use ligand ASL:", args.lig_asl

    num_frames_processed = 0
    for frame_index in xrange(start_frame, end_frame, args.step_size):
        print "Reading frame %i..." % frame_index
        sys.stdout.flush()
        st = csim.getFrameStructure(frame_index)

        try:
            del st.property['s_chorus_trajectory_file']
        except KeyError:
            pass

        st.write('complex.mae')

        if water_atoms is None:
            water_atoms = analyze.evaluate_asl(st, '(res.ptype "CLR ") OR (res.ptype "OLA ") OR (res.ptype "OLC ") OR (res.ptype "OLB ") OR (res.ptype "POPC") OR (res.ptype "SPC ") OR (res.ptype " NA ") OR (res.ptype "CL  ")')
        st.deleteAtoms(water_atoms)

        # Find the ligand:
        if args.lig_asl is None:
            args.lig_asl = find_ligand(st)
            print "Will use ligand ASL:", args.lig_asl

        writer.append(st)
        num_frames_processed += 1
    writer.close()

    if args.res_file:
        print 'Running Residue Scanning on %i complexes' % num_frames_processed
        complexes_file = run_residue_scanning(complexes_file,
                                              args.res_file,
                                              args.lig_asl,
                                              args.jobname)

    num_prime_sts = structure.count_structures(complexes_file)
    print 'Passing %i structures to Prime' % num_prime_sts
    print "Running Prime MMGBSA job..."

    # Run a multi-job MMGBSA calculation:
    prime_jobname = args.jobname + '-prime'
    cmd = [PRIME_MMGBSA, '-job_type', 'REAL_MIN', complexes_file,
           '-jobname', prime_jobname, '-ligand', args.lig_asl]

    # Determine what -HOST value the user has specified:
    host_str = jobcontrol.host_list_to_str(jobcontrol.get_command_line_host_list())
    if not args.NJOBS:
        # Use the number of CPUs:
        args.NJOBS = jobcontrol.calculate_njobs()
    cmd += ['-NJOBS', str(args.NJOBS), '-HOST', host_str, '-OVERWRITE']

    print "Running:", subprocess.list2cmdline(cmd)

    job = jobcontrol.launch_job(cmd)
    print "  Prime log file: %s.log" % prime_jobname
    sys.stdout.flush()
    job.wait()

    if not job.succeeded():
        print "ERROR: Prime MMGBSA job failed. JobID: %s" % job.job_id
        sys.exit(1)

    print "Parsing Prime MMGBSA output file..."
    sys.stdout.flush()
    prime_out_file = prime_jobname + "-out.maegz"
    prime_csv_file = prime_jobname + "-out.csv"
    print 'Got %i structures from Prime' % structure.count_structures(prime_out_file)

    props_to_report = MMGBSA_PROPS.items()
    if args.res_file:
        props_to_report += RESIDUE_SCANNING_PROPS.items()
        mut_strings = []

    values_by_prop = defaultdict(list)
    for complex_st in structure.StructureReader(prime_out_file):
        for _, dataname in props_to_report:
            value = complex_st.property[dataname]
            values_by_prop[dataname].append(value)
        if args.res_file:
            mut_strings.append(complex_st.property[MUT_PROP])

    # Report average and standard deviation for each property:
    for name, dataname in props_to_report:
        values = values_by_prop[dataname]
        ave, std = calc_ave_and_std(values)
        print "%s Average: %.4f" % (name, ave)
        print "%s Standard Deviation: %.2f" % (name, std)
        print "%s Range: %.4f to %.4f" % (name, min(values), max(values))
        print ""

    print "Number of frames processed: %i" % num_frames_processed
    print "Number of output structures: %i" % len(values)
    print "Output structure file:", prime_out_file
    print "Prime CSV file:", prime_csv_file
    if args.res_file:
        # Generate the CSV file containing per-mutation result data.
        muts_csv_file = args.jobname + '-mutations.csv'
        muts_csv_fh = open(muts_csv_file, "wb")
        muts_csv_writer = csv.writer(muts_csv_fh)
        header = ["Mutation"]
        for name, _ in props_to_report:
            header += ["Ave_" + name, "StdDev_" + name, "Range_" + name]
        muts_csv_writer.writerow(header)

        # Group outputs by mutation:
        indices_by_mut = OrderedDict()  # key: mut_str; value: index to
        # values list
        for i, mut_str in enumerate(mut_strings):
            try:
                indices_by_mut[mut_str].append(i)
            except KeyError:
                indices_by_mut[mut_str] = [i]

        # Calculate property averages per mutation:
        for mut_str, mut_indices in indices_by_mut.iteritems():
            row = [mut_str]
            for name, dataname in props_to_report:
                all_values = values_by_prop[dataname]
                mut_values = [all_values[index] for index in mut_indices]
                ave, std = calc_ave_and_std(mut_values)
                range_str = "%.2f to %.2f" % (min(mut_values), max(mut_values))
                row += [round(ave, 4), round(std, 2), range_str]
            muts_csv_writer.writerow(row)
        muts_csv_fh.close()
        print "Mutation results CSV file:", muts_csv_file
    sys.stdout.flush()
