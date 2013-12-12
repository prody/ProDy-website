# -*- coding: utf-8 -*-

import os
import sys
import glob
import time

try:
    exec(open('manual/conf.py').read())
except IOError:
    exec(open('../../manual/conf.py').read())

sys.path.append(os.path.abspath('manual/sphinxext'))
              #'sphinxcontrib.googleanalytics',
              #'ipython_console_highlighting',
              #'ipython_directive']

exclude_patterns.append('ProDy')
exclude_patterns.append('tutorials/template')
exclude_patterns.extend(glob.glob('tutorials/**/acknowledgments.rst'))
templates_path = ['_template']
source_suffix = '.rst'
master_doc = 'contents'

# not needed when building the full website
intersphinx_mapping.pop('prodywebsite')

project = u'ProDy'
copyright = u'2010-2014, University of Pittsburgh'

# -- Options for HTML output ---------------------------------------------------

html_favicon = 'manual/_static/favicon.ico'

html_additional_pages = {'index': 'index.html'}

# -- Options for LaTeX output --------------------------------------------------
latex_logo = 'manual/_static/logo.png'

lines = (line for line in rst_epilog.split('\n')
         if ('_Tut' not in line and '_NMW' not in line))
rst_epilog = '\n'.join(lines)