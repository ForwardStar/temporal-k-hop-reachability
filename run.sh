echo "Compiling..."
make
echo "Running..."
if [ -n "$3" ]; then
    ./main graph.txt query.txt output.txt $1 $2 $3
elif [ -n "$2" ]; then
    ./main graph.txt query.txt output.txt $1 $2
elif [ -n "$1" ]; then
    ./main graph.txt query.txt output.txt $1
else
    ./main graph.txt query.txt output.txt AdvancedTwoHopIndex BFS
fi