#export ROOTCOREGRID=1 && $1/RootCore/scripts/grid_compile_config.sh $1 && source $1/RootCore/scripts/setup.sh && $ROOTCOREDIR/scripts/compile.sh
export ROOTCOREGRID=1 && $1/RootCore/scripts/grid_compile_config.sh $1 && source $1/RootCore/scripts/setup.sh && $ROOTCOREDIR/scripts/grid_spec_compile.sh
