"""=========================================Pipeline single cell test data generation============================================Overview=========This pipeline implements Splatter and Minnow to create simulated 10X single cell read fastq files.Requires:* a SCE (single cell experiment object, saved in RDS format)- to estimate parameters from    * Saved as format <details>_sce.rds* yml file for number of genes, number of cells, number of count matrices producedPipeline output================* Numerous simulated count matrices* Numerous simulated fastq filesCode====="""from ruffus import *#import sys#import osimport cgatcore.pipeline as Pimport cgatcore.experiment as E# load options from the config filePARAMS = P.get_parameters(    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],     "../pipeline.yml",     "pipeline.yml"])INPUT_SCE = ["*.sce.rds"]@follows(mkdir("out_dir"))@split(INPUT_SCE,       "out.dir/*/quants_mat.csv")def splatter_generate(infile, outfiles):    """    Runs Rscript to generate count matrices, gene names and cell barcode csv files    """
    # Add optparse for command line arguments    statement = """ Rscript splatter_simulate_channel_patients.R     """

    job_memory = "100G"

    P.run(statement)@follows(splatter_generate)
@transform(   etc.  ,

)
def minnow_run(infiles, outfiles):
    """
    Runs minnow in splatter mode.
    """ 


    statement = """ minnow --splatter ....
    """

    job_memory = "100G"

    P.run(statement)


@follows(splatter_generate, minnow_run)def full():    pass

def main(argv=None):    if argv is None:        argv = sys.argv    P.main(argv)if __name__ == "__main__":    sys.exit(P.main(sys.argv))


