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

config=`pwd`"/mof-config/run-y3vo2-mof.yaml"
tiles="test-tiles.yaml"

megamix-meds                       \
    --skip-errors                  \
    --noblind                      \
    $config                          \
    $cmd                           \
    $tiles
