CXX = g++
CXXFLAGS = -Wall -std=c++17
LDFLAGS =
SOURCE_DIR = src
INCLUDE_DIR = includes
TEST_DIR = tests

BUILD_DIR = build
BIN_DIR = bin
ALL_BUILD_DIR = $(BUILD_DIR) $(BIN_DIR)

TEST_BUILD_DIR = $(TEST_DIR)/build
TEST_BIN_DIR = $(TEST_DIR)/bin
ALL_TEST_BUILD_DIR = $(TEST_BUILD_DIR) $(TEST_BIN_DIR)

all: SPH_2D 

SPH_2D: $(BIN_DIR)/SPH_2D

$(BIN_DIR)/SPH_2D: $(BUILD_DIR)/SPH_2D.o $(BUILD_DIR)/SPH_Snippet.o $(BUILD_DIR)/file_writer.o
	$(CXX) -o $@ $^

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.cpp $(INCLUDE_DIR)/*.h | directories
	$(CXX) $(CPPFLAGS) -o $@ -c $< $(CXXFLAGS) -I$(INCLUDE_DIR)

clean:
	rm -f $(BUILD_DIR)/* $(BIN_DIR)/* *.vtp

.PHONY: SPH_2D all clean

TESTS = test_SPH_2D test_file_writer

runtests: cleantest ${TESTS}
	@python3 run_tests.py

tests: ${TESTS}

test_SPH_2D: $(TEST_BIN_DIR)/test_SPH_2D

test_file_writer: $(TEST_BIN_DIR)/test_file_writer

$(TEST_BIN_DIR)/test_SPH_2D: $(TEST_BUILD_DIR)/test_SPH_2D.o $(BUILD_DIR)/SPH_2D.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS)

$(TEST_BIN_DIR)/test_file_writer: $(TEST_BUILD_DIR)/test_file_writer.o $(BUILD_DIR)/file_writer.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS)

$(TEST_BUILD_DIR)/%.o: $(TEST_DIR)/%.cpp $(INCLUDE_DIR)/*.h | test_directories
	$(CXX) -o $@ -c $< $(CXXFLAGS) $(CPPFLAGS) -I$(INCLUDE_DIR)

cleantest:
	rm -f $(TEST_BUILD_DIR)/* $(TEST_BIN_DIR)/* $(TEST_DIR)/*.vtp

.PHONY: tests ${TESTS} cleantests runtests


directories:
	@mkdir -p $(ALL_BUILD_DIR)

test_directories:
	@mkdir -p $(ALL_TEST_BUILD_DIR)

.PHONY: directories test_directories
