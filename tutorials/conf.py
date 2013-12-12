exec(open('../../conf.py').read())
sys.path.append(os.path.abspath('../../manual/sphinxext'))

# needed when building PDF files for tutorials separately
intersphinx_mapping['prodywebsite'] = ('http://prody.csb.pitt.edu', None)

master_doc = 'index'

version = release = tutorial_version or version

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
