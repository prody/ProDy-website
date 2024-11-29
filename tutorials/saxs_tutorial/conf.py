tutorial_title = u'SAXS Analysis'
tutorial_author = u'Mustafa Tekpinar'
tutorial_logo = u''           # default is ProDy logo
tutorial_version = u''        # default is latest ProDy version

# keep the following part as is
# comment out this -> we don't compile sax
try:
    execfile('../conf.py')
except NameError:
    exec(open('../conf.py').read())
