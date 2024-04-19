import os
import fypp
import argparse
from joblib import Parallel, delayed

C_PREPROCESSED = ( 

)

def pre_process_toml(kargs):
    """Pre process mannifest
    fpm.toml mannifest file can be edited to account for custom keywords

    Parameters
    ----------
    kargs : str
        collection of command line arguments.
    """
    from tomlkit import table, dumps
    data = table()
    data.add("name", "fsparse")
    data.add("version", str(kargs.vmajor)+
                    "."+str(kargs.vminor)+
                    "."+str(kargs.vpatch) )
    data.add("license", "MIT")
    data.add("author", "José Alves")
    data.add("maintainer", "jose.alves@transvalor.com")
    data.add("copyright", "Copyright 2023, Jose Alves")

    dev_dependencies = table()
    dev_dependencies.add("test-drive", {"git" : "https://github.com/fortran-lang/test-drive", 
                                        "tag" : "v0.4.0"})
    data.add("dev-dependencies", dev_dependencies)

    preprocess = table()
    preprocess.add("cpp", {} )
    preprocess['cpp'].add("suffixes", [".f90"] )
    # preprocess['cpp'].add("macros", [] )
    data.add("preprocess", preprocess)

    with open("fpm.toml", "w") as f:
        f.write(dumps(data))

def pre_process_fypp(kargs):
    """Pre process fypp files

    Parameters
    ----------
    kargs : str
        collection of command line arguments.
    """
    kwd = []
    if kargs.with_qp:
        kwd.append("-DWITH_QP=True")
    if kargs.with_xqp:
        kwd.append("-DWITH_XQP=True")
    
    optparser = fypp.get_option_parser()
    options, leftover = optparser.parse_args(args=kwd)
    options.includes = ['include'] #folder with include files for fypp
    # options.line_numbering = True
    tool = fypp.Fypp(options)

    # Define the folders to search for *.fypp files
    folders = ['src', 'test']
    # Process all folders
    fypp_files = [os.path.join(root, file) for folder in folders
              for root, _, files in os.walk(folder)
              for file in files if file.endswith(".fypp")]
    
    def process_f(file):
        source_file = file
        basename = os.path.splitext(source_file)[0]
        sfx = 'f90' # if os.path.basename(basename) not in C_PREPROCESSED else 'F90'
        target_file = basename + '.' + sfx
        tool.process_file(source_file, target_file)
    
    Parallel(n_jobs=kargs.njob)(delayed(process_f)(f) for f in fypp_files)
    
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Preprocess FSPARSE source files.')
    # fypp arguments
    parser.add_argument("--vmajor", type=int, default=0, help="Project Version Major")
    parser.add_argument("--vminor", type=int, default=1, help="Project Version Minor")
    parser.add_argument("--vpatch", type=int, default=1, help="Project Version Patch")

    parser.add_argument("--njob", type=int, default=4, help="number of parallel jobs")
    parser.add_argument("--with_qp", action="store_true", help="Include WITH_QP in the command")
    parser.add_argument("--with_xqp", action="store_true", help="Include WITH_XQP in the command")
    parser.add_argument("--with_hp", action="store_true", help="Include WITH_HP in the command")
    
    args = parser.parse_args()

    pre_process_toml(args)
    pre_process_fypp(args)