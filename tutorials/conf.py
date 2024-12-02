import sys
import os
intersphinx_mapping = {
    'python': ('http://docs.python.org/', None),
    'numpy': ('http://docs.scipy.org/doc/numpy/', None),
    'scipy': ('http://docs.scipy.org/doc/scipy/reference/', None),
    'matplotlib': ('http://matplotlib.sourceforge.net/', None),
    'prodywebsite': ('http://www.bahargroup.org/prody/', None),
}
RTD = os.environ.get('READTHEDOCS', None) == 'True'
intersphinx_mapping.pop('prodywebsite')
if RTD:
    html_theme = 'default'
else:
    templates_path = ['_theme']
    html_theme = '_theme'
    html_theme_path = ['.']
    html_static_path = ['_static']
html_theme_options = {}

extensions = [ 'sphinx.ext.todo',
               'sphinx.ext.autodoc',
               'sphinx.ext.doctest',
               'sphinx.ext.coverage',
               'sphinx.ext.extlinks',
               'sphinx.ext.graphviz',
               'sphinx.ext.ifconfig',
               'sphinx.ext.viewcode',
               'sphinx.ext.intersphinx',
               'sphinx.ext.inheritance_diagram',
               'matplotlib.sphinxext.mathmpl',
               'matplotlib.sphinxext.only_directives',
               'IPython.sphinxext.ipython_directive' ,
               'IPython.sphinxext.ipython_console_highlighting',
               ]

#exec(open('../conf.py').read())
sys.path.append(os.path.abspath('../manual/sphinxext'))

# needed when building PDF files for tutorials separately
intersphinx_mapping['prodywebsite'] = ('http://www.bahargroup.org/prody', None)
master_doc = 'index'

#version = release = tutorial_version #or version
latex_documents = [
    ('index',
     os.path.basename(os.getcwd()) + '.tex',
     tutorial_title,
     tutorial_author,
     'manual'),
]

latex_logo = tutorial_logo or '_static/logo.png'
latex_show_urls = 'footnote'

html_additional_pages = {}
html_domain_indices = False
html_use_index = False

latex_domain_indices = True

latex_appendices = ['acknowledgments']
