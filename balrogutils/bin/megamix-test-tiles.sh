#   Possible commands are:
#        setup - setup jobs
#        setup-nbrs - setup jobs for neighbors/fof finding
#        collate - combine job outputs into a single file
#        verify - verify that all job outputs are present and OK
#        clean - clean all outputs from a run
#        archive - run after collate to delete intermediate files and tar logs

if [[ $# -lt 1 ]]; then
    echo "usage: $(basename $0) cmd"
    exit 1
fi

cmd=$1

if [[ $cmd == "setup-missing" ]]; then
    cmd="setup"
    extra="--missing"
    echo "setup for missing"
else
    extra=""
fi

# you could make your own config and use that here
config=`pwd`"/sof-config/run-y3v02-sof.yaml"
#config="/data/des71.a/data/kuropat/blank_test/balrog-base/desmeds-config/meds-y3v02.yaml"
tiles="test-tiles.yaml"

megamix-meds                       \
    --skip-errors                  \
    --noblind                      \
    $config                        \
    $cmd                           \
    $tiles
