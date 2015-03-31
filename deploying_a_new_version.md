# Why this way? #

Currently there is no way of building directly on the archive machine itself (no gcc), so the source must first be compiled on ltdevsrv. The resulting binaries are then transferred to lt-archive via scp.

# Method #

## Compiling the source ##

  * ssh on to ltdevsrv

> ` [eng@rmb-tower]$ ssh eng@ltdevsrv `

  * checkout the latest version from svn

> ` [eng@ltdevsrv ~/rmb]$ svn checkout http://frodo-l2-pipeline.googlecode.com/svn/trunk/ frodo-l2-pipeline-read-only `

  * alter the L2\_BASE\_DIR variable in /scripts correspondingly

> ` [eng@ltdevsrv scripts]$ vi L2_setup `

  * source L2\_setup

> ` [eng@ltdevsrv scripts]$ source L2_setup `

  * alter Makefile in src/ and uncomment ltdevsrv associated LIBS and INCLUDES variables, commenting any others

> ` [eng@ltdevsrv src]$ vi Makefile `

  * issue make all command from src/

> ` [eng@ltdevsrv src]$ make all `

  * copy the resulting binaries in bin/ to lt-archive

> ` [eng@ltdevsrv bin]$ scp * data@lt-archive:"{PATH}" `

  * remove files from ltdevsrv

> ` [eng@ltdevsrv ~/rmb]$ rm -rf frodo-l2-pipeline-read-only/ `

## Deploying on lt-archive ##

  * ssh on to lt-archive

> ` [eng@rmb-tower]$ ssh data@lt-archive `

  * checkout the latest version from svn

> ` [data@lt-archive ~/frodo_build]$ svn checkout http://frodo-l2-pipeline.googlecode.com/svn/ frodo_l2_pipeline `

  * move the binaries compiled from ltdevsrv to bin/ directory

> ` [data@lt-archive ~/frodo_build]$ mv * frodo_l2_pipeline/trunk/bin/ `

  * alter the L2\_BASE\_DIR, L2\_SUCCESS\_LOG\_PATH, L2\_FAIL\_LOG\_PATH variables in /scripts correspondingly

> ` [eng@lt-archive scripts]$ vi L2_setup `

  * as root, suffix old frodo directory with bak and move new version in its place

> ` [root@lt-archive bin]# su `

> ` [root@lt-archive bin]# mv frodo_l2_pipeline/ frodo_l2_pipeline_bak `

> ` [root@lt-archive bin]# mv {PATH}/frodo_l2_pipeline . `

  * move previous config directories (reference\_arcs, lookup\_tables) to new directories

> ` [root@lt-archive config]# mv {PATH}/frodo_l2_pipeline_bak/trunk/config/reference_arcs/ . `

> ` [root@lt-archive config]# mv {PATH}/frodo_l2_pipeline_bak/trunk/config/lookup_tables/ . `