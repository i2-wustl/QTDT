#
# QTDT Makefile -- Compiles and installs QTDT and accessory applications
# (c) 1999-2005 Goncalo Abecasis
#

# version information
VERSION = 2.6.1
PSVERSION = 0.6.12

# default installation directory
INSTALLDIR=/usr/local/bin

# default C++ compiler and flags
CXX=g++
CFLAGS=-O -I./libsrc -D__ZLIB_AVAILABLE__

# executable file names and locations
BINDIR = executables
QTDT = $(BINDIR)/qtdt
PRELUDE = $(BINDIR)/prelude
FINALE = $(BINDIR)/finale
PEDSTATS = $(BINDIR)/pedstats
EXECUTABLES = $(QTDT) $(PRELUDE) $(FINALE) $(PEDSTATS)

# QTDT File Set
QTDTSRC = qtdt/Qtdt.cpp qtdt/QtdtUi.cpp qtdt/QtdtDescribe.cpp \
 qtdt/QtdtLinear.cpp qtdt/QtdtVar.cpp qtdt/QtdtOptions.cpp \
 qtdt/QtdtBroker.cpp qtdt/QtdtParameters.cpp qtdt/QtdtResult.cpp
QTDTHDR = $(QTDTSRC:.cpp=.h)
QTDTOBJ = $(QTDTSRC:.cpp=.o)

# Utility Library File Set
LIBFILE = libsrc/lib-goncalo.a
LIBMAIN = libsrc/BasicHash libsrc/Error libsrc/FortranFormat \
 libsrc/GenotypeLists libsrc/IBD \
 libsrc/InputFile libsrc/IntArray libsrc/Hash \
 libsrc/LongArray libsrc/Kinship libsrc/KinshipX  libsrc/MapFunction \
 libsrc/MathCholesky libsrc/MathDeriv libsrc/MathFloatVector \
 libsrc/MathGenMin libsrc/MathGold libsrc/MathMatrix libsrc/MathStats \
 libsrc/MathNormal libsrc/MathSVD libsrc/MathVector \
 libsrc/MemoryInfo libsrc/MiniDeflate \
 libsrc/Parameters libsrc/Pedigree libsrc/PedigreeAlleleFreq \
 libsrc/PedigreeDescription libsrc/PedigreeFamily libsrc/PedigreeGlobals \
 libsrc/PedigreePerson libsrc/QuickIndex libsrc/Random libsrc/Sort \
 libsrc/StringArray libsrc/StringBasics libsrc/StringMap \
 libsrc/StringHash libsrc/TraitTransformations
LIBPED = libsrc/PedigreeLoader libsrc/PedigreeTwin libsrc/PedigreeTrim
LIBSRC = $(LIBMAIN:=.cpp) $(LIBPED:=.cpp)
LIBHDR = $(LIBMAIN:=.h) libsrc/Constant.h libsrc/Genetics.h \
 libsrc/MathConstant.h libsrc/PedigreeAlleles.h libsrc/LongInt.h
LIBOBJ = $(LIBSRC:.cpp=.o)

# Other Files
QTDTEXAMPLES = examples

# private parameters
FETCHDIR=$(HOME)/code
DISTRIBDIR=$(HOME)/code/distrib/qtdt-$(VERSION)

# helpful screen listing available options
help : 
	@echo "QTDT Source Distribution"
	@echo " "
	@echo "This Makefile will compile and install qtdt on your system"
	@echo " "
	@echo "Type...           To..."
	@echo "make help         Display this help screen"
	@echo "make all          Compile qtdt, prelude, finale and pedstats"
	@echo "make install      Install binaries in $(INSTALLDIR)"
	@echo "make install INSTALLDIR=directory_for_binaries"
	@echo "                  Install binaries in directory_for_binaries"
	@echo "make clean        Delete temporary files"

# make everything
all : $(EXECUTABLES)

$(EXECUTABLES) : $(BINDIR)

$(BINDIR) :
	mkdir $(BINDIR)

# dependencies for executables
$(QTDT) : $(LIBFILE) $(QTDTOBJ)
	$(CXX) $(CFLAGS) -o $@ $(QTDTOBJ) $(LIBFILE) -lm -lz

$(PRELUDE) : $(LIBFILE) extras/prelude.cpp $(LIBHDR)
	$(CXX) $(CFLAGS) -o $@ extras/prelude.cpp $(LIBFILE) -lz 

$(FINALE) : $(LIBFILE) extras/finale.cpp $(LIBHDR)
	$(CXX) $(CFLAGS) -o $@ extras/finale.cpp $(LIBFILE) -lz

$(PEDSTATS) : pedstats-$(PSVERSION).tar.gz
	gunzip -c pedstats-$(PSVERSION).tar.gz | tar -xf -
	cd pedstats-$(PSVERSION) ; $(MAKE) executables/pedstats
	cp pedstats-$(PSVERSION)/executables/pedstats executables
	rm -rf pedstats-$(PSVERSION)

$(LIBFILE) : $(LIBOBJ) $(LIBHDR)
	ar -cr $@ $(LIBOBJ)
	ranlib $@

$(QTDTOBJ) : $(QTDTHDR) $(LIBHDR)

$(LIBOBJ) : $(LIBHDR)

clean :
	-rm -f */*.a */*.o $(EXECUTABLES) 

install : all $(INSTALLDIR)
	@echo Installing to directory $(INSTALLDIR)
	@echo To select a different directory, run
	@echo " "
	@echo make install INSTALLDIR=your_preferred_dir
	@echo " "
	@cp $(EXECUTABLES) $(INSTALLDIR)

$(INSTALLDIR) :
	@echo Creating directory $(INSTALLDIR)
	@echo " "
	@mkdir $(INSTALLDIR)

new-version :
	mkdir -p $(DISTRIBDIR) $(DISTRIBDIR)/extras 
	mkdir -p $(DISTRIBDIR)/qtdt $(DISTRIBDIR)/libsrc
	cp ChangeLog LICENSE.twister README $(DISTRIBDIR)
	cp Makefile $(DISTRIBDIR)
	cp -R $(QTDTEXAMPLES) $(DISTRIBDIR)

fetch : 
	cd $(FETCHDIR) ; cp $(QTDTSRC) $(QTDTHDR) $(DISTRIBDIR)/qtdt
	cd $(FETCHDIR) ; cp $(LIBSRC) $(LIBHDR) $(DISTRIBDIR)/libsrc
	cp $(FETCHDIR)/prelude/prelude.cpp $(DISTRIBDIR)/extras
	cp $(FETCHDIR)/prelude/finale.cpp $(DISTRIBDIR)/extras
	cd $(DISTRIBDIR) ; wget -r -nd -N http://www.sph.umich.edu/csg/abecasis/pedstats/download/pedstats-$(PSVERSION).tar.gz
	cd $(DISTRIBDIR) ; csh ../stamp QTDT

.c.o :
	$(CXX) $(CFLAGS) -o $@ -c $*.c

.cpp.o : 
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp -DVERSION=\"$(VERSION)\"

archive : clean
	mkdir -p qtdt-$(VERSION)
	cp -R ChangeLog LICENSE.twister README Makefile qtdt-$(VERSION)
	cp -R libsrc qtdt extras $(QTDTEXAMPLES) qtdt-$(VERSION)
	cp -R pedstats-$(PSVERSION).tar.gz qtdt-$(VERSION)
	tar -cvf qtdt-$(VERSION).tar qtdt-$(VERSION) 
	gzip -f --best qtdt-$(VERSION).tar
	rm -rf qtdt-$(VERSION)
	tar -cvf examples.tar $(QTDTEXAMPLES)
	gzip -f --best examples.tar

distrib : $(EXECUTABLES)
	mkdir -p qtdt-$(VERSION)
	cp -R ChangeLog LICENSE.twister README $(EXECUTABLES) $(QTDTEXAMPLES) qtdt-$(VERSION)
	tar -cvf `uname`-qtdt.tar qtdt-$(VERSION)
	gzip -f `uname`-qtdt.tar
	rm -rf qtdt-$(VERSION)

windowszip : $(EXECUTABLES)
	mkdir -p qtdt-$(VERSION)
	cp -R ChangeLog LICENSE.twister README $(EXECUTABLES) $(QTDTEXAMPLES) qtdt-$(VERSION)
	echo '@ECHO OFF' > qtdt-$(VERSION)/STARTHERE.bat
	echo 'CMD /K "SET PATH=%CD%;%PATH%"' > qtdt-$(VERSION)/STARTHERE.bat
	zip -r Windows-qtdt.zip qtdt-$(VERSION)
	rm -rf qtdt-$(VERSION)

.SUFFIXES : .cpp .c .o $(SUFFIXES)

