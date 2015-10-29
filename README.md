frodo-l2-pipeline
=================

# Installation on lt-qc

* grab a copy of the dist ($git clone https://github.com/LivTel/frodo-l2-pipeline.git)
* edit ROOT/scripts/L2\_setup and change L2\_BASE\_DIR pointer to lt-qc
* set up the environment (source ROOT/scripts/L2\_setup)
* edit ROOT/src/Makefile to point $LIBS and $INCLUDES to lt-qc
* make binaries and library (ROOT/src/make all)





