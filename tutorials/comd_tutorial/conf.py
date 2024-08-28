tutorial_title = u'coMD Tutorial'
tutorial_author = u'Cihan Kaya'
tutorial_logo = u''           # default is ProDy logo
tutorial_version = u''        # default is latest ProDy version
extensions = ['IPython.sphinxext.ipython_console_highlighting',
              'IPython.sphinxext.ipython_directive']
# keep the following part as is
try:
    execfile('../conf.py')
except NameError:
    exec(open('../conf.py').read())
