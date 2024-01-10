#include "Graph.h"
#include "popl.h"
#include "Utility.h"

#include <string>
#include <csignal>

using namespace std;
using namespace popl;
void print_usage()
{
    printf("Example usage: ./kDefective -g path_to_graph -k 3 -o\n");
}
FILE *f;
static void SIGINT_exit(int signum)
{
    printf("Time Out!\n");
    fprintf(f, "NA\tNA\tNA\tNA\tNA\tNA\n ");
    fclose(f);
    exit(1);
}

int main(int argc, char *argv[])
{
    bool output = false;

    OptionParser op("Allowed options");
    auto help_option = op.add<Switch>("h", "help", "\'produce help message\'");
    auto graph_option = op.add<Value<string>>("g", "graph", "\'path to input graph file\'");
    auto k_option = op.add<Value<int>>("k", "k", "\'the value of k for k-plex\'");
    op.add<Switch>("o", "output", "\'write the kplex to ./kdefective.txt\'", &output);

    op.parse(argc, argv);

    if (help_option->is_set() || argc <= 1)
    {
        cout << op << endl;
        if (argc <= 1)
        {
            print_usage();
            return 0;
        }
    }
    if (!graph_option->is_set())
    {
        printf("!!! Path to input graph file is not provided! Exit !!!\n");
        return 0;
    }
    if (!k_option->is_set())
    {
        printf("!!! k is not provided! Exit !!!\n");
        return 0;
    }
    f = fopen("result.txt", "a+");
    fprintf(f, "%s\t%d\t", graph_option->value().c_str(), k_option->value());
    Graph *graph = new Graph(graph_option->value().c_str(), k_option->value());
    signal(SIGTERM, SIGINT_exit);
    graph->read_graph();
    graph->k_defective_exact(f);
    if (!graph->verify_result())
    {
        printf("!!!verify failed\n");
        return 0;
    };

    graph->output_result();
    if (output)
    {
        graph->output_one_kdefective();
    }
    fclose(f);
    return 0;
}