import re

SMK_EXTENSIONS= ['.smk', '.snake', '.snakefile']

def read_snakefile(snakefile):
    with open(snakefile) as fin:
        smk= fin.read()
    return smk

def get_docstrings(snakestring):
    """Parse the snakefile and return a list of dictionaries. Each dictionary 
    has key inidcating the type of element being documented and value the 
    markdown string.
    """
    docs= re.findall('"""<.+?"""', snakestring, flags= re.DOTALL)
    docs= [re.sub('^"""<\s+', '', x) for x in docs]
    docs= [re.sub('"""$', '', x) for x in docs]
    return docs

def write_readme(docstrings, readme_file):
    """Write the docstrings to the readme file
    """
    with open(readme_file, 'w') as fout:
        for x in docstrings:
            fout.write(x + '\n')

def compose_readme_fname(snakefile, smk_exts= SMK_EXTENSIONS, suffix= '.README.md'):
    """Return a filename for the README based on the given 
    snakefile template.
    """
    outname= None
    for x in smk_exts:
        if snakefile.endswith(x):
            outname= re.sub(x + '$', '', snakefile) + suffix
            break
    if outname is None:
        outname= snakefile + suffix
    return outname
    
def make_readme(snakefile, readme_file= None):
    """Main entry point to write out the readme file based the input snakemake file
    """
    docs= get_docstrings(read_snakefile(snakefile))
    if readme_file is None:
        readme_file= compose_readme_fname(snakefile)
    write_readme(docs, readme_file)
