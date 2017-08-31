NVCC = nvcc
SRC = src
EXEC_DIR = bin
OBJ_DIR = obj
EXEC = $(EXEC_DIR)/swift
TEST_EXEC = $(EXEC_DIR)/test
OUTPUT_TEST_EXEC = $(EXEC_DIR)/testOutput
CUDA_ARCH = sm_20
PROFILE = -pg -arch=$(CUDA_ARCH)
DEBUG = -g -G -gencode arch=compute_20,code=sm_20 -gencode arch=compute_10,code=sm_10 
NVCCFLAGS = -lpthread -O3 -arch=$(CUDA_ARCH)
TEST = test
CHECK_LIBS = -lcheck
CHECK_LIBS_PATH = -L/usr/lib
CHECK_FLAGS = $(CHECK_LIBS_PATH) $(CHECK_LIBS)
MAIN_SOURCE_FILE = $(SRC)/main.cu
SOURCE_FILES = $(SRC)/preprocess.cu \
				$(SRC)/input.cu \
				$(SRC)/search.cu \
				$(SRC)/search2.cu \
				$(SRC)/query.cu \
				$(SRC)/reference.cu \
				$(SRC)/refPosList.cu \
				$(SRC)/hitList.cu \
				$(SRC)/align.cu \
				$(SRC)/memory.cu \
				$(SRC)/output.cu \
				$(SRC)/array.cu \
				$(SRC)/smithWaterman.cu \
				$(SRC)/refPosMap.cu \
				$(SRC)/refNameMap.cu \
				$(SRC)/refMap.cu \
				$(SRC)/lookupTable.cu \
				$(SRC)/lookupTable2.cu \
				$(SRC)/lookupTable3.cu \
				$(SRC)/lookupTable4.cu \
				$(SRC)/lookupTable5.cu \
				$(SRC)/lookupTable6.cu \
				$(SRC)/lookupTable7.cu \
				$(SRC)/lookupTable8.cu \
				$(SRC)/lookupTable9.cu \
				$(SRC)/lookupTable10.cu \
				$(SRC)/lookupTable11.cu \
				$(SRC)/lookupTable13.cu \
				$(SRC)/mapHits.cu \
				$(SRC)/mapHits2.cu \
				$(SRC)/mapHits3.cu \
				$(SRC)/mapHits4.cu \
				$(SRC)/mapHits5.cu \
				$(SRC)/mapHits6.cu

TEST_SOURCE_FILES = $(TEST)/testHitList.cu \
				$(TEST)/testRefPosList.cu \
				$(TEST)/testReference.cu \
				$(TEST)/testArray.cu \
				$(TEST)/testInput.cu \
				$(TEST)/testQuery.cu \
				$(TEST)/testSearch.cu \
				$(TEST)/testAlign.cu \
				$(TEST)/testSmithWaterman.cu \
				$(TEST)/testRefPosMap.cu \
				$(TEST)/testRefNameMap.cu \
				$(TEST)/testRefMap.cu \
				$(TEST)/testLookupTable.cu \
				$(TEST)/testLookupTable5.cu \
				$(TEST)/testLookupTable6.cu \
				$(TEST)/testLookupTable7.cu \
				$(TEST)/testMapHits.cu \
				$(TEST)/testMapHits2.cu \
				$(TEST)/testMapHits3.cu \
				$(TEST)/testMapHits4.cu \
				$(TEST)/testMapHits5.cu \
				$(TEST)/testMapHits6.cu \
				$(TEST)/testMain.cu
				
OUTPUT_TEST_FILES = $(TEST)/testAlignOutput.cu
				
MAIN_OBJ_FILE = $(OBJ_DIR)/main.o
				
OBJ_FILES = $(OBJ_DIR)/preprocess.o \
				$(OBJ_DIR)/input.o \
				$(OBJ_DIR)/search.o \
				$(OBJ_DIR)/search2.o \
				$(OBJ_DIR)/query.o \
				$(OBJ_DIR)/reference.o \
				$(OBJ_DIR)/refPosList.o \
				$(OBJ_DIR)/hitList.o \
				$(OBJ_DIR)/align.o \
				$(OBJ_DIR)/memory.o \
				$(OBJ_DIR)/output.o \
				$(OBJ_DIR)/array.o \
				$(OBJ_DIR)/smithWaterman.o \
				$(OBJ_DIR)/refPosMap.o \
				$(OBJ_DIR)/refNameMap.o \
				$(OBJ_DIR)/refMap.o \
				$(OBJ_DIR)/lookupTable.o \
				$(OBJ_DIR)/lookupTable2.o \
				$(OBJ_DIR)/lookupTable3.o \
				$(OBJ_DIR)/lookupTable4.o \
				$(OBJ_DIR)/lookupTable5.o \
				$(OBJ_DIR)/lookupTable6.o \
				$(OBJ_DIR)/lookupTable7.o \
				$(OBJ_DIR)/lookupTable8.o \
				$(OBJ_DIR)/lookupTable9.o \
				$(OBJ_DIR)/lookupTable10.o \
				$(OBJ_DIR)/lookupTable11.o \
				$(OBJ_DIR)/lookupTable13.o \
				$(OBJ_DIR)/mapHits.o \
				$(OBJ_DIR)/mapHits2.o \
				$(OBJ_DIR)/mapHits3.o \
				$(OBJ_DIR)/mapHits4.o \
				$(OBJ_DIR)/mapHits5.o \
				$(OBJ_DIR)/mapHits6.o
				
TEST_OBJ_FILES = $(OBJ_DIR)/testHitList.o \
				$(OBJ_DIR)/testRefPosList.o \
				$(OBJ_DIR)/testReference.o \
				$(OBJ_DIR)/testArray.o \
				$(OBJ_DIR)/testInput.o \
				$(OBJ_DIR)/testQuery.o \
				$(OBJ_DIR)/testSearch.o \
				$(OBJ_DIR)/testAlign.o \
				$(OBJ_DIR)/testSmithWaterman.o \
				$(OBJ_DIR)/testRefPosMap.o \
				$(OBJ_DIR)/testRefNameMap.o \
				$(OBJ_DIR)/testRefMap.o \
				$(OBJ_DIR)/testLookupTable.o \
				$(OBJ_DIR)/testLookupTable5.o \
				$(OBJ_DIR)/testLookupTable6.o \
				$(OBJ_DIR)/testLookupTable7.o \
				$(OBJ_DIR)/testMapHits.o \
				$(OBJ_DIR)/testMapHits2.o \
				$(OBJ_DIR)/testMapHits3.o \
				$(OBJ_DIR)/testMapHits4.o \
				$(OBJ_DIR)/testMapHits5.o \
				$(OBJ_DIR)/testMapHits6.o \
				$(OBJ_DIR)/testMain.o

				
all: $(MAIN_OBJ_FILE) $(OBJ_FILES)
	$(NVCC) $(MAIN_OBJ_FILE) $(OBJ_FILES) -o $(EXEC)
	
check: $(OBJ_FILES) $(TEST_OBJ_FILES)
	$(NVCC) $(CHECK_FLAGS) $(OBJ_FILES) $(TEST_OBJ_FILES) -o $(TEST_EXEC)
	
testOutput: 
	$(NVCC) $(NVCCFLAGS) $(OUTPUT_TEST_FILES) -o $(OUTPUT_TEST_EXEC)
	
profile:
	$(NVCC) $(PROFILE) $(MAIN_SOURCE_FILE) $(SOURCE_FILES) -o $(EXEC)
	
debug:
	$(NVCC) $(DEBUG) $(MAIN_SOURCE_FILE) $(SOURCE_FILES) -o $(EXEC)
	
emu:
	$(NVCC) -deviceemu $(DEBUG) $(MAIN_SOURCE_FILE) $(SOURCE_FILES) \
	-o $(EXEC_DIR)/emu
	
clean:
	rm -fr $(OBJ_DIR)/*.o $(EXEC) $(TEST_EXEC)
	
$(OBJ_DIR)/main.o: $(SRC)/main.cu $(SRC)/preprocess.h $(SRC)/input.h \
$(SRC)/common.h $(SRC)/search.h $(SRC)/search2.h $(SRC)/align.h $(SRC)/refMap.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/main.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/preprocess.o: $(SRC)/preprocess.cu $(SRC)/preprocess.h \
$(SRC)/common.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/preprocess.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/input.o: $(SRC)/input.cu $(SRC)/common.h $(SRC)/input.h 
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/input.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/search.o: $(SRC)/search.cu $(SRC)/search.h $(SRC)/common.h \
$(SRC)/input.h $(SRC)/preprocess.h $(SRC)/reference.h $(SRC)/lookupTable.h \
$(SRC)/lookupTable2.h $(SRC)/lookupTable3.h $(SRC)/lookupTable4.h \
$(SRC)/lookupTable5.h $(SRC)/lookupTable6.h $(SRC)/lookupTable7.h \
$(SRC)/lookupTable8.h $(SRC)/lookupTable9.h $(SRC)/lookupTable10.h \
$(SRC)/lookupTable11.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/search.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/search2.o: $(SRC)/search2.cu $(SRC)/search2.h $(SRC)/common.h \
$(SRC)/input.h $(SRC)/preprocess.h 
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/search2.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/query.o: $(SRC)/query.cu $(SRC)/query.h $(SRC)/common.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/query.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/reference.o: $(SRC)/reference.cu $(SRC)/reference.h \
$(SRC)/common.h $(SRC)/preprocess.h $(SRC)/array.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/reference.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/refPosList.o: $(SRC)/refPosList.cu $(SRC)/refPosList.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/refPosList.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/hitList.o: $(SRC)/hitList.cu $(SRC)/hitList.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/hitList.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/align.o: $(SRC)/align.cu $(SRC)/align.h $(SRC)/common.h \
$(SRC)/smithWaterman.h $(SRC)/memory.h $(SRC)/output.h $(SRC)/input.h \
$(SRC)/query.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/align.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/memory.o: $(SRC)/memory.cu $(SRC)/memory.h $(SRC)/common.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/memory.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/output.o: $(SRC)/output.cu $(SRC)/output.h $(SRC)/common.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/output.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/array.o: $(SRC)/array.cu $(SRC)/array.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/array.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/smithWaterman.o: $(SRC)/smithWaterman.cu $(SRC)/smithWaterman.h \
$(SRC)/common.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/smithWaterman.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/refPosMap.o: $(SRC)/refPosMap.cu $(SRC)/refPosMap.h \
$(SRC)/common.h $(SRC)/refPosList.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/refPosMap.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/refNameMap.o: $(SRC)/refNameMap.cu $(SRC)/refNameMap.h \
$(SRC)/common.h 
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/refNameMap.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/refMap.o: $(SRC)/refMap.cu $(SRC)/refMap.h $(SRC)/common.h 
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/refMap.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/lookupTable.o: $(SRC)/lookupTable.cu $(SRC)/lookupTable.h \
$(SRC)/common.h $(SRC)/preprocess.h $(SRC)/array.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/lookupTable.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/lookupTable2.o: $(SRC)/lookupTable2.cu $(SRC)/lookupTable2.h \
$(SRC)/common.h $(SRC)/preprocess.h $(SRC)/array.h $(SRC)/refMap.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/lookupTable2.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/lookupTable3.o: $(SRC)/lookupTable3.cu $(SRC)/lookupTable3.h \
$(SRC)/common.h $(SRC)/preprocess.h $(SRC)/array.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/lookupTable3.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/lookupTable4.o: $(SRC)/lookupTable4.cu $(SRC)/lookupTable4.h \
$(SRC)/common.h $(SRC)/preprocess.h $(SRC)/mapHits.h $(SRC)/mapHits2.h \
$(SRC)/mapHits3.h $(SRC)/mapHits4.h $(SRC)/mapHits5.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/lookupTable4.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/lookupTable5.o: $(SRC)/lookupTable5.cu $(SRC)/lookupTable5.h \
$(SRC)/common.h $(SRC)/preprocess.h $(SRC)/mapHits6.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/lookupTable5.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/lookupTable6.o: $(SRC)/lookupTable6.cu $(SRC)/lookupTable6.h \
$(SRC)/common.h $(SRC)/preprocess.h $(SRC)/mapHits6.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/lookupTable6.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/lookupTable7.o: $(SRC)/lookupTable7.cu $(SRC)/lookupTable7.h \
$(SRC)/common.h $(SRC)/preprocess.h $(SRC)/mapHits6.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/lookupTable7.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/lookupTable8.o: $(SRC)/lookupTable8.cu $(SRC)/lookupTable8.h \
$(SRC)/common.h $(SRC)/preprocess.h $(SRC)/mapHits6.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/lookupTable8.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/lookupTable9.o: $(SRC)/lookupTable9.cu $(SRC)/lookupTable9.h \
$(SRC)/common.h $(SRC)/preprocess.h $(SRC)/mapHits6.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/lookupTable9.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/lookupTable10.o: $(SRC)/lookupTable10.cu $(SRC)/lookupTable10.h \
$(SRC)/common.h $(SRC)/preprocess.h $(SRC)/mapHits6.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/lookupTable10.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/lookupTable11.o: $(SRC)/lookupTable11.cu $(SRC)/lookupTable11.h \
$(SRC)/common.h $(SRC)/preprocess.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/lookupTable11.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/lookupTable13.o: $(SRC)/lookupTable13.cu $(SRC)/lookupTable13.h \
$(SRC)/common.h $(SRC)/preprocess.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/lookupTable13.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/mapHits.o: $(SRC)/mapHits.cu $(SRC)/mapHits.h $(SRC)/common.h \
$(SRC)/array.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/mapHits.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/mapHits2.o: $(SRC)/mapHits2.cu $(SRC)/mapHits2.h $(SRC)/common.h \
$(SRC)/array.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/mapHits2.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/mapHits3.o: $(SRC)/mapHits3.cu $(SRC)/mapHits3.h $(SRC)/common.h \
$(SRC)/array.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/mapHits3.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/mapHits4.o: $(SRC)/mapHits4.cu $(SRC)/mapHits4.h $(SRC)/common.h \
$(SRC)/array.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/mapHits4.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/mapHits5.o: $(SRC)/mapHits5.cu $(SRC)/mapHits5.h $(SRC)/common.h \
$(SRC)/array.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/mapHits5.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/mapHits6.o: $(SRC)/mapHits6.cu $(SRC)/mapHits6.h $(SRC)/common.h \
$(SRC)/array.h
	$(NVCC) $(NVCCFLAGS) -c $(SRC)/mapHits6.cu -odir $(OBJ_DIR)
    
$(OBJ_DIR)/testLookupTable.o: $(TEST)/testLookupTable.cu $(SRC)/lookupTable.h \
$(TEST)/testLookupTable.h $(SRC)/lookupTable2.h $(SRC)/lookupTable3.h \
$(SRC)/refMap.h $(SRC)/preprocess.h 
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testLookupTable.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testLookupTable5.o: $(TEST)/testLookupTable5.cu $(SRC)/lookupTable5.h \
$(TEST)/testLookupTable5.h $(SRC)/refMap.h $(SRC)/preprocess.h 
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testLookupTable5.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testLookupTable6.o: $(TEST)/testLookupTable6.cu $(SRC)/lookupTable6.h \
$(TEST)/testLookupTable6.h $(SRC)/refMap.h $(SRC)/preprocess.h 
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testLookupTable6.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testLookupTable7.o: $(TEST)/testLookupTable7.cu $(SRC)/lookupTable7.h \
$(TEST)/testLookupTable7.h $(SRC)/refMap.h $(SRC)/preprocess.h 
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testLookupTable7.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testMapHits.o: $(TEST)/testMapHits.cu $(SRC)/mapHits.h \
$(TEST)/testMapHits.h  
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testMapHits.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testMapHits2.o: $(TEST)/testMapHits2.cu $(SRC)/mapHits2.h \
$(TEST)/testMapHits2.h  
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testMapHits2.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testMapHits3.o: $(TEST)/testMapHits3.cu $(SRC)/mapHits3.h \
$(TEST)/testMapHits3.h  
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testMapHits3.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testMapHits4.o: $(TEST)/testMapHits4.cu $(SRC)/mapHits4.h \
$(TEST)/testMapHits4.h  
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testMapHits4.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testMapHits5.o: $(TEST)/testMapHits5.cu $(SRC)/mapHits5.h \
$(TEST)/testMapHits5.h  
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testMapHits5.cu -odir $(OBJ_DIR)

$(OBJ_DIR)/testMapHits6.o: $(TEST)/testMapHits6.cu $(SRC)/mapHits6.h \
$(TEST)/testMapHits6.h  
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testMapHits6.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testRefMap.o: $(TEST)/testRefMap.cu $(SRC)/refMap.h \
$(TEST)/testRefMap.h 
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testRefMap.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testRefNameMap.o: $(TEST)/testRefNameMap.cu $(SRC)/refNameMap.h \
$(TEST)/testRefNameMap.h 
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testRefNameMap.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testRefPosMap.o: $(TEST)/testRefPosMap.cu $(SRC)/refPosMap.h \
$(TEST)/testRefPosMap.h $(SRC)/reference.h
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testRefPosMap.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testHitList.o: $(TEST)/testHitList.cu $(SRC)/hitList.h \
$(TEST)/testHitList.h
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testHitList.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testRefPosList.o: $(TEST)/testRefPosList.cu \
$(TEST)/testRefPosList.h $(SRC)/refPosList.h
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testRefPosList.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testReference.o: $(TEST)/testReference.cu $(TEST)/testReference.h \
$(SRC)/reference.h $(SRC)/preprocess.h
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testReference.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testArray.o: $(TEST)/testArray.cu $(TEST)/testArray.h $(SRC)/array.h
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testArray.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testInput.o: $(TEST)/testInput.cu $(TEST)/testInput.h \
$(SRC)/input.h $(SRC)/common.h
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testInput.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testQuery.o: $(TEST)/testQuery.cu $(TEST)/testQuery.h \
$(SRC)/query.h $(SRC)/common.h
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testQuery.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testSearch.o: $(TEST)/testSearch.cu $(TEST)/testSearch.h \
$(SRC)/search.h $(SRC)/refMap.h $(SRC)/preprocess.h
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testSearch.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testAlign.o: $(TEST)/testAlign.cu $(TEST)/testAlign.h \
$(SRC)/align.h $(SRC)/common.h
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testAlign.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testSmithWaterman.o: $(TEST)/testSmithWaterman.cu \
$(TEST)/testSmithWaterman.h $(SRC)/smithWaterman.h $(SRC)/common.h
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testSmithWaterman.cu -odir $(OBJ_DIR)
	
$(OBJ_DIR)/testMain.o: $(TEST)/testMain.cu $(TEST)/testHitList.h \
$(TEST)/testRefPosList.h $(TEST)/testReference.h $(TEST)/testArray.h \
$(TEST)/testInput.h $(TEST)/testQuery.h $(TEST)/testSearch.h \
$(TEST)/testAlign.h $(TEST)/testSmithWaterman.h $(TEST)/testRefPosMap.h \
$(TEST)/testRefNameMap.h $(TEST)/testRefMap.h $(TEST)/testLookupTable.h \
$(TEST)/testLookupTable5.h $(TEST)/testLookupTable6.h $(TEST)/testLookupTable6.h \
$(TEST)/testMapHits.h $(TEST)/testMapHits2.h $(TEST)/testMapHits4.h \
$(TEST)/testMapHits5.h $(TEST)/testMapHits6.h
	$(NVCC) $(CHECK_FLAGS) -c $(TEST)/testMain.cu -odir $(OBJ_DIR)
