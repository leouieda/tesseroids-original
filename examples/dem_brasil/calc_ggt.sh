#/bin/bash

export TESSBINPATH="../../bin/"

${TESSBINPATH}tessgrd -r-60/-45/-30/-15 -b100/100 -z250e03 | ${TESSBINPATH}tessgxx $1 -a | \
${TESSBINPATH}tessgxy $1 -a | ${TESSBINPATH}tessgxz $1 -a | ${TESSBINPATH}tessgyy $1 -a | \
${TESSBINPATH}tessgyz $1 -a | ${TESSBINPATH}tessgzz $1 -a > $2