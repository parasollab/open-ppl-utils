TREE_TYPE = TWO_LEVEL_TREE
# TREE_TYPE = THREE_LEVEL_TREE
# TREE_TYPE = ONE_LEVEL_TREE

CC = g++
IDIR = INCLUDE/LKHInclude
ODIR = OBJ
CFLAGS = -O3 -Wall -I$(IDIR) -D$(TREE_TYPE) -g

_DEPS = Delaunay.h GainType.h Genetic.h GeoConversion.h Hashing.h      \
        Heap.h LKH.h Segment.h Sequence.h

DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = Activate.o AddCandidate.o AddExtraCandidates.o                  \
       AddTourCandidates.o AdjustCandidateSet.o                        \
       AllocateStructures.o Ascent.o                                   \
       Best2OptMove.o Best3OptMove.o Best4OptMove.o Best5OptMove.o     \
       BestKOptMove.o Between.o Between_SL.o Between_SSL.o             \
       BridgeGain.o BuildKDTree.o C.o CandidateReport.o                \
       ChooseInitialTour.o Connect.o CreateCandidateSet.o              \
       CreateDelaunayCandidateSet.o CreateQuadrantCandidateSet.o       \
       Delaunay.o Distance.o Distance_SPECIAL.o eprintf.o ERXT.o       \
       Excludable.o Exclude.o FindTour.o FixedOrCommonCandidates.o     \
       Flip.o Flip_SL.o Flip_SSL.o Forbidden.o FreeStructures.o        \
       fscanint.o Gain23.o GenerateCandidates.o Genetic.o              \
       GeoConversion.o GetTime.o GreedyTour.o Hashing.o Heap.o         \
       IsBackboneCandidate.o IsCandidate.o IsCommonEdge.o              \
       IsPossibleCandidate.o KSwapKick.o LinKernighan.o                \
       Make2OptMove.o Make3OptMove.o Make4OptMove.o Make5OptMove.o     \
       MakeKOptMove.o MergeTourWithBestTour.o MergeWithTour.o          \
       Minimum1TreeCost.o MinimumSpanningTree.o NormalizeNodeList.o    \
       NormalizeSegmentList.o OrderCandidateSet.o PatchCycles.o        \
       printff.o PrintParameters.o                                     \
       Random.o ReadCandidates.o ReadLine.o ReadParameters.o           \
       ReadPenalties.o ReadProblem.o RecordBestTour.o                  \
       RecordBetterTour.o RemoveFirstActive.o                          \
       ResetCandidateSet.o RestoreTour.o SegmentSize.o Sequence.o      \
       SFCTour.o SolveCompressedSubproblem.o                           \
       SolveDelaunaySubproblems.o SolveKarpSubproblems.o               \
       SolveKCenterSubproblems.o SolveKMeansSubproblems.o              \
       SolveRoheSubproblems.o SolveSFCSubproblems.o SolveSubproblem.o  \
       SolveSubproblemBorderProblems.o SolveTourSegmentSubproblems.o   \
       Statistics.o StoreTour.o SymmetrizeCandidateSet.o               \
       TrimCandidateSet.o WriteCandidates.o WritePenalties.o           \
       WriteTour.o LKH.o

OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)


standalone:
	$(MAKE) LKH

LKH: $(OBJ) $(DEPS)
	$(CC) -o ../LKH $(OBJ) $(CFLAGS) -lm

NAME   := LKH
STATIC := lib$(NAME).a

.DEFAULT_GOAL := $(STATIC)

$(STATIC): $(OBJ) | $(DEPS)
	@echo Linking static library...
	@ar rcs ../$@ $^

.PHONY: clean
clean:
	/bin/rm -f $(ODIR)/*.o ../LKH *~ ._* $(IDIR)/*~ $(IDIR)/._* ../$(STATIC)

