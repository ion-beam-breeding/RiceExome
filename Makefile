# common settings
include Makefile.in

## project-wide settings ##
BINDIR = bin
CXXFLAGS +=
LDLIBS +=
INCLUDES +=

## multi-file targets ##
TARGETS = adjust_target_location genotype_filter vcf2xls pindel_vcf_filter merge_vc \
bt_coverage_filter

## adjust_target_location ##
SRCS_adjust_target_location = adjust_target_location.cpp histd.cpp seq.cpp
OBJS_adjust_target_location = $(SRCS_adjust_target_location:.cpp=.o)
CFLAGS_adjust_target_location =
LDLIBS_adjust_target_location =

## genotype_filter ##
SRCS_genotype_filter = genotype_filter.cpp histd.cpp
OBJS_genotype_filter = $(SRCS_genotype_filter:.cpp=.o)
CFLAGS_genotype_filter =
LDLIBS_genotype_filter =

## vcf2xls ##
SRCS_vcf2xls = vcf2xls.cpp histd.cpp
OBJS_vcf2xls = $(SRCS_vcf2xls:.cpp=.o)
CFLAGS_vcf2xls =
LDLIBS_vcf2xls =

## pindel_vcf_filter ##
SRCS_pindel_vcf_filter = pindel_vcf_filter.cpp histd.cpp
OBJS_pindel_vcf_filter = $(SRCS_pindel_vcf_filter:.cpp=.o)
CFLAGS_pindel_vcf_filter =
LDLIBS_pindel_vcf_filter =

## merge_vc ##
SRCS_merge_vc = merge_vc.cpp histd.cpp
OBJS_merge_vc = $(SRCS_merge_vc:.cpp=.o)
CFLAGS_merge_vc =
LDLIBS_merge_vc =

## bt_coverage_filter ##
SRCS_bt_coverage_filter = bt_coverage_filter.cpp histd.cpp
OBJS_bt_coverage_filter = $(SRCS_bt_coverage_filter:.cpp=.o)
CFLAGS_bt_coverage_filter =
LDLIBS_bt_coverage_filter =



## common targets ##
all: $(BINDIR) $(TARGETS) $(MISCPROGS)
clean:
	$(RM) $(BINDIR)/*
	$(RM) *.o
$(BINDIR):
	$(MKDIR) -p $(BINDIR)
$(OBJS_$@): Makefile
.cpp.o:
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_$@) $(INCLUDES) $(INCLUDES_$@) -c $< -o $@
.c.o:
	$(CC) $(CFLAGS) $(CFLAGS_$@) $(INCLUDES) $(INCLUDES_$@) -c $< -o $@
.cpp:
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_$@) $(INCLUDES) $(INCLUDES_$@) $< -o $(BINDIR)/$@ $(LDLIBS) $(LDLIBS_$@)

## adjust_target_location ##
adjust_target_location: $(BINDIR) $(OBJS_adjust_target_location)
	$(LD) $(CXXFLAGS) $(CXXFLAGS_$@) -o $(BINDIR)/$@ $(OBJS_$@) $(LDLIBS) $(LDLIBS_$@)
adjust_target_location_clean:
	$(RM) $(BINDIR)/adjust_target_location

## genotype_filter ##
genotype_filter: $(BINDIR) $(OBJS_genotype_filter)
	$(LD) $(CXXFLAGS) $(CXXFLAGS_$@) -o $(BINDIR)/$@ $(OBJS_$@) $(LDLIBS) $(LDLIBS_$@)
genotype_filter_clean:
	$(RM) $(BINDIR)/genotype_filter

## vcf2xls ##
vcf2xls: $(BINDIR) $(OBJS_vcf2xls)
	$(LD) $(CXXFLAGS) $(CXXFLAGS_$@) -o $(BINDIR)/$@ $(OBJS_$@) $(LDLIBS) $(LDLIBS_$@)
vcf2xls_clean:
	$(RM) $(BINDIR)/vcf2xls

## pindel_vcf_filter ##
pindel_vcf_filter: $(BINDIR) $(OBJS_pindel_vcf_filter)
	$(LD) $(CXXFLAGS) $(CXXFLAGS_$@) -o $(BINDIR)/$@ $(OBJS_$@) $(LDLIBS) $(LDLIBS_$@)
pindel_vcf_filter_clean:
	$(RM) $(BINDIR)/pindel_vcf_filter

## merge_vc ##
merge_vc: $(BINDIR) $(OBJS_merge_vc)
	$(LD) $(CXXFLAGS) $(CXXFLAGS_$@) -o $(BINDIR)/$@ $(OBJS_$@) $(LDLIBS) $(LDLIBS_$@)
merge_vc_clean:
	$(RM) $(BINDIR)/merge_vc

## bt_coverage_filter ##
bt_coverage_filter: $(BINDIR) $(OBJS_bt_coverage_filter)
	$(LD) $(CXXFLAGS) $(CXXFLAGS_$@) -o $(BINDIR)/$@ $(OBJS_$@) $(LDLIBS) $(LDLIBS_$@)
bt_coverage_filter_clean:
	$(RM) $(BINDIR)/bt_coverage_filter


## dependency check ##
.KEEP_STATE:
.KEEP_STATE_FILE:.make.state.GNU-x86-Linux
