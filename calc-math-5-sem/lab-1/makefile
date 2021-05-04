# ====================
# COMPILER OPTIONS

CXX = g++

CXXFLAGS = -Wall -g

PY = python3

# ====================
# VARIABLES

# directories
BUILD_DIR  	 = build
OUTPUT_DIR 	 = output
SOURCE_DIR 	 = src

CALC_DIR_SRC = ${SOURCE_DIR}/calculator
DRAW_DIR_SRC = ${SOURCE_DIR}/drawer

CALC_DIR_OBJ = ${BUILD_DIR}/calculator
DRAW_DIR_OBJ = ${BUILD_DIR}/drawer

CALC_DIR_OUT = ${OUTPUT_DIR}/calculator
DRAW_DIR_OUT = ${OUTPUT_DIR}/drawer

#libraries
CALC_OBJ_FILES = ${patsubst ${CALC_DIR_SRC}/%.cpp, ${CALC_DIR_OBJ}/%.o, ${wildcard ${CALC_DIR_SRC}/*.cpp}}

DRAW_SRC_FILES = ${wildcard ${DRAW_DIR_SRC}/*.py}

# ====================
# RULES

all: clear calc draw

# drawer
draw: ${DRAW_SRC_FILES} | ${DRAW_DIR_OUT}
	${PY} $^

# calculator
calc: calc.out | ${CALC_DIR_OUT}
	./${BUILD_DIR}/calc.out

calc.out: ${CALC_OBJ_FILES}
	${CXX} ${CXXFLAGS} -o ${BUILD_DIR}/$@ $^

${CALC_DIR_OBJ}/%.o: ${CALC_DIR_SRC}/%.cpp | ${CALC_DIR_OBJ} 
	${CXX} ${CXXFLAGS} -c $< -o $@

clear:
	clear

clean:
	rm -rf build/ output/ calc.out

# directory gen
${BUILD_DIR}:
	mkdir ${BUILD_DIR}

${CALC_DIR_OBJ}:| ${BUILD_DIR}
	mkdir ${BUILD_DIR}/calculator

${DRAW_DIR_OBJ}: | ${BUILD_DIR}
	mkdir ${BUILD_DIR}/drawer

${OUTPUT_DIR}:
	mkdir ${OUTPUT_DIR}

${CALC_DIR_OUT}: | ${OUTPUT_DIR}
	mkdir ${OUTPUT_DIR}/calculator

${DRAW_DIR_OUT}: | ${OUTPUT_DIR}
	mkdir ${OUTPUT_DIR}/drawer




