# ====================
# COMPILER OPTIONS

CXX = g++

CXXFLAGS = -Wall -g

PY = python3

# ====================
# VARIABLES

# directories
BUILD_DIR  = build
OUTPUT_DIR = output
SOURCE_DIR = src

CALC_LIB_DIR_SRC = ${SOURCE_DIR}/calculator/lib
CALC_DIR_SRC     = ${SOURCE_DIR}/calculator
DRAW_DIR_SRC     = ${SOURCE_DIR}/drawer

CALC_LIB_DIR_OBJ = ${BUILD_DIR}/calculator/lib
CALC_DIR_OBJ     = ${BUILD_DIR}/calculator

CALC_DIR_OUT = ${OUTPUT_DIR}/calculator
DRAW_DIR_OUT = ${OUTPUT_DIR}/drawer

# object files
CALC_LIB_OBJ_FILES = ${patsubst ${CALC_LIB_DIR_SRC}/%.cpp, ${CALC_LIB_DIR_OBJ}/%.o, ${wildcard ${CALC_LIB_DIR_SRC}/*.cpp}}
CALC_LIB_HDR_FILES = ${wildcard ${CALC_LIB_DIR_SRC}/*.h}

CALC_OBJ_FILES = ${patsubst ${CALC_DIR_SRC}/%.cpp, ${CALC_DIR_OBJ}/%.o, ${wildcard ${CALC_DIR_SRC}/*.cpp}}
CALC_HDR_FILES = ${wildcard ${CALC_DIR_SRC}/*.h}

# source files
DRAW_SRC_FILES = ${wildcard ${DRAW_DIR_SRC}/*.py}

# ====================
# RULES

all: CALC

# drawer
DRAW: ${DRAW_SRC_FILES} | ${DRAW_DIR_OUT}
	${PY} $^

# calculator
CALC: CALC.out | ${CALC_DIR_OUT}

CALC.out: ${CALC_LIB_OBJ_FILES} ${CALC_OBJ_FILES} ${CALC_LIB_HDR_FILES}
	${CXX} ${CXXFLAGS} -o ${BUILD_DIR}/$@ $^

${CALC_DIR_OBJ}/%.o: ${CALC_DIR_SRC}/%.cpp | ${CALC_DIR_OBJ} 
	${CXX} ${CXXFLAGS} -c $< -o $@

# matrix
${CALC_LIB_DIR_OBJ}/%.o: ${CALC_LIB_DIR_SRC}/%.cpp | ${CALC_LIB_DIR_OBJ} 
	${CXX} ${CXXFLAGS} -c $< -o $@

clear:
	clear

clean:
	rm -rf build/ output/

# directory gen
${BUILD_DIR}:
	mkdir ${BUILD_DIR}

${CALC_LIB_DIR_OBJ}: | ${CALC_DIR_OBJ}
	mkdir ${BUILD_DIR}/calculator/lib

${CALC_DIR_OBJ}: | ${BUILD_DIR}
	mkdir ${BUILD_DIR}/calculator

${OUTPUT_DIR}:
	mkdir ${OUTPUT_DIR}

${CALC_DIR_OUT}: | ${OUTPUT_DIR}
	mkdir ${OUTPUT_DIR}/calculator

${DRAW_DIR_OUT}: | ${OUTPUT_DIR}
	mkdir ${OUTPUT_DIR}/drawer




