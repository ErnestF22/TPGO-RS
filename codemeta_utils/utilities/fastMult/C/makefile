OBJDIR=build
BINDIR=bin
MEXDIR=mex

CFLAGS:=-I. -Wall -Wextra -O3 -fPIC
CC:=gcc

MEXBIN=$(shell $(MEXUTILITYDIR)/detectMexBin)
MEXUTILITYDIR:=../../mexUtility
MEXOPT:=-I$(MEXUTILITYDIR)
MEX:=$(MEXBIN) $(MEXOPT)
MEXEXT:=$(shell $(MEXUTILITYDIR)/detectMexExt)


.PHONY: all allmex allLibraries allTests clean fastMult arrayUtility 

OBJS=fastMult arrayUtility
MEXS= mexFastMult3x3MatMat mexFastMult3x3HatMat mexFastMult3x3MatHat mexFastMult3x3HatNormSqMat mexFastMult3x3MatHatNormSq mexFastMult3x3MatMatMat mexFastMult3x3MatIdxTranspMatMatIdx

#add prefix and suffix to object and mex files
OBJS:=$(addsuffix .o, $(OBJS))
OBJS:=$(addprefix $(OBJDIR)/, $(OBJS))
MEXS:=$(addsuffix .$(MEXEXT), $(MEXS))
MEXS:=$(addprefix $(MEXDIR)/, $(MEXS))


all: arrayUtility fastMult allmex

clean:
	-rm -rf $(BINDIR)/*.dSYM 
	-rm $(OBJDIR)/*.o $(BINDIR)/* $(MEXDIR)/*.$(MEXEXT)

#libraries
allLibraries: $(OBJS)

fastMult: $(OBJDIR)/fastMult.o $(BINDIR)/fastMult_test

arrayUtility: $(OBJDIR)/arrayUtility.o $(BINDIR)/arrayUtility_test


#generic object files
$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

#test executables
allTests: $(BINDIR)/fastMult_test $(BINDIR)/arrayUtility_test

#generic executables
$(BINDIR)/% : %.c  $(OBJS)
	$(CC) $(CFLAGS) $*.c $(OBJDIR)/*.o -o $@

#MEX FUNCTIONS
#allmex: $(MEXDIR)/mexFastMult3x3MatMat.$(MEXEXT) $(MEXDIR)/mexFastMult3x3HatMat.$(MEXEXT) $(MEXDIR)/mexFastMult3x3MatHat.$(MEXEXT) $(MEXDIR)/mexFastMult3x3HatNormSqMat.$(MEXEXT) $(MEXDIR)/mexFastMult3x3MatHatNormSq.$(MEXEXT) $(MEXDIR)/mexFastMult3x3MatMatMat.$(MEXEXT)
allmex: $(MEXS)

$(MEXDIR)/%.$(MEXEXT): $(MEXDIR)/%.c $(OBJDIR)/fastMult.o
	$(MEX) $< $(OBJDIR)/*.o -outdir $(MEXDIR)
