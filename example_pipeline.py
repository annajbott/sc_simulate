"""=========================================
    # Add optparse for command line arguments

    job_memory = "100G"

    P.run(statement)
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


@follows(splatter_generate, minnow_run)

def main(argv=None):

