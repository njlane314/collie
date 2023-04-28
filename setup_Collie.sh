source /pc2014-data2/pguzowski/snsw_3.0/setup.sh
source /pc2014-data2/pguzowski/snsw_3.0/CadfaelBrew/Cellar/root6/6.08.06/libexec/thisroot.sh

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
COLLIE="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
#COLLIE=`dirname $BASH_SOURCE`
export LD_LIBRARY_PATH="$COLLIE/lib:$COLLIE/CLHEP/lib:$LD_LIBRARY_PATH"
export PATH="$COLLIE/examples:$PATH"
#ln -sf /dev/null examples/fort.99
