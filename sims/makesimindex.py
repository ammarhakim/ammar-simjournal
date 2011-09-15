from mako.template import Template
myTemplate = Template(filename='simindex.template')
fp = open("simindex.rst", "w")
fp.writelines( myTemplate.render() )
