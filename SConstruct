# Build the Tesseroids programs
import sys

# get the mode flag from the command line
mode = ARGUMENTS.get('mode', 'default')
if not (mode in ['default', 'check', 'bin32', 'win32']):
   print "Error: unknown mode '%s'" % (mode)
   Exit(1)
print '**** Compiling in ' + mode + ' mode...'

if sys.platform == 'win32':
    env = Environment(
        CPPPATH='src/lib')
elif mode == 'check':
    env = Environment(
        CFLAGS='-ansi -pedantic-errors -Wall -ggdb',
        LIBS=['m'],
        CPPPATH='src/lib')
elif mode == 'win32':
    env = Environment(
        CFLAGS='-O3',
        LIBS=['m'],
        CPPPATH='src/lib')
    env.Tool('crossmingw', toolpath=['scons-tools'])
elif mode == 'bin32':
    env = Environment(
        CFLAGS='-O3 -m32',
        LINKFLAGS='-m32',
        LIBS=['m'],
        CPPPATH='src/lib')
else:
    env = Environment(
        CFLAGS='-O3',
        LIBS=['m'],
        CPPPATH='src/lib')

# Build the tessg* programs
tesssrc = Split("""
    src/lib/logger.c
    src/lib/version.c
    src/lib/grav_tess.c
    src/lib/glq.c
    src/lib/constants.c
    src/lib/geometry.c
    src/lib/parsers.c
    src/lib/tessg_main.c
    """)
fields = ['pot', 'gx', 'gy', 'gz', 'gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz']
for f in fields:
    sources = ['src/tess%s.c' % (f)] + tesssrc
    env.Program('bin/tess%s' % (f), source=sources)

# Build the prismg* programs
tesssrc = Split("""
    src/lib/logger.c
    src/lib/version.c
    src/lib/constants.c
    src/lib/geometry.c
    src/lib/parsers.c
    src/lib/prismg_main.c
    src/lib/grav_prism.c
    """)
fields = ['pot', 'gx', 'gy', 'gz', 'gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz']
for f in fields:
    sources = ['src/prism%s.c' % (f)] + tesssrc
    env.Program('bin/prism%s' % (f), source=sources)

# Build prismpots, prismgs, and prismggts
env.Program('bin/prismpots', source=Split("""
    src/prismpots.c
    src/lib/grav_prism_sph.c
    src/lib/grav_prism.c
    src/lib/logger.c
    src/lib/version.c
    src/lib/constants.c
    src/lib/geometry.c
    src/lib/parsers.c
    """))
env.Program('bin/prismgs', source=Split("""
    src/prismgs.c
    src/lib/grav_prism_sph.c
    src/lib/grav_prism.c
    src/lib/logger.c
    src/lib/version.c
    src/lib/constants.c
    src/lib/geometry.c
    src/lib/parsers.c
    """))
env.Program('bin/prismggts', source=Split("""
    src/prismggts.c
    src/lib/grav_prism_sph.c
    src/lib/grav_prism.c
    src/lib/logger.c
    src/lib/version.c
    src/lib/constants.c
    src/lib/geometry.c
    src/lib/parsers.c
    """))

# Build tess2prism
env.Program('bin/tess2prism', source=Split("""
    src/tess2prism.c
    src/lib/logger.c
    src/lib/version.c
    src/lib/constants.c
    src/lib/geometry.c
    src/lib/parsers.c
    """))
# Build tessdefaults
env.Program('bin/tessdefaults', source=Split("""
    src/tessdefaults.c
    src/lib/logger.c
    src/lib/version.c
    src/lib/constants.c
    src/lib/glq.c
    src/lib/geometry.c
    """))
# Build tessgrd
env.Program('bin/tessgrd', source=Split("""
    src/tessgrd.c
    src/lib/logger.c
    src/lib/version.c
    src/lib/parsers.c
    src/lib/constants.c
    """))
# Build tessmass
env.Program('bin/tessmass', source=Split("""
    src/tessmass.c
    src/lib/logger.c
    src/lib/version.c
    src/lib/parsers.c
    src/lib/geometry.c
    src/lib/constants.c
    """))
# Build tessmodgen
env.Program('bin/tessmodgen', source=Split("""
    src/tessmodgen.c
    src/lib/logger.c
    src/lib/version.c
    src/lib/parsers.c
    src/lib/geometry.c
    src/lib/constants.c
    """))
# Build tesslayers
env.Program('bin/tesslayers', source=Split("""
    src/tesslayers.c
    src/lib/logger.c
    src/lib/version.c
    src/lib/parsers.c
    src/lib/geometry.c
    src/lib/constants.c
    """))

# Build the test runner
sources = ['test/test_all.c']
sources.extend(Glob("src/lib/*.c"))
tesstest = env.Program('tesstest', source=sources)

# Clean exe files
Clean('.', Glob('bin/*.exe'))
Clean('.', Glob('*.exe'))
