/*

This is a Graph class implementation.

Version 2.0

Date of Last modified: 25/11/2021

Last modified by: Kirill Arestov


*/


#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <ctime>
#include <iterator>


using namespace std;


class Graph {
private:

	// Class that stores a vertex info. Internal class of class Graph
	class Vertex {
	private:

	public:
		const int INIT_COST = 1000000;

		int Cost, Name;
		bool is_destination;				// Var for detection whether a vertex is a destination for Dijkstra's algorithm

		// Constructor by default
		Vertex()
		{
			Name = 0;
			Cost = INIT_COST;
			is_destination = 0;
		}

		int GetInitCost() { return INIT_COST; }

	};


	// VARIABLES AND FUNCTIONS OF CLASS GRAPH (private)

	// General variables
	const int
		MIN_NUM_OF_VERTEX = 2,
		MAX_NUM_OF_VERTEX = 410;

	//Graph Gra;

	vector<vector<int>> MainGraph;								// Vector NxN for storing a graph in matrix orientation
	Vertex* Vertexes = new Vertex[MAX_NUM_OF_VERTEX + 1];		// Dynamic array of Vertex objects. Stores all the vertexes in a graph. Initializing with maximum available size for limiting overheads	
	int
		num_of_vertexes = 0,									// Num. of vertexes in a graph
		num_of_edges = 0;										// Num. of edges in a graph

	// Variables for Dijkstra's algorithm
	vector <int> Path, ShortestPath;							// Vectors for testing and accumulating final closed sets
	int path_counter = 0;

	// Variables for MST construction algorithms
	vector <Graph> SimpleTrees;																						// Vector for storing simple trees (all edges) for Kruskal's algorithm
	int num_of_simple_trees = 0;


	// GENERAL FUNCTIONS
	//
	// Creates zero graph with 2 vertexes. Uses for some constructors
	void create_zero_graph()
	{
		num_of_vertexes = 2;

		MainGraph.resize(num_of_vertexes + 1, vector<int>(num_of_vertexes + 1));								// Creates graph in NxN matrix representation

		for (int i = 1; i < num_of_vertexes; i++) {
			MainGraph[0][i] = MainGraph[i][0] = Vertexes[i].Name = i;											// Generate the names of each vertex (from 1 to num_of_vertexes)

			for (int j = i + 1; j <= num_of_vertexes; j++) {													// Generate 0 values for all edges by default
				MainGraph[i][j] = MainGraph[j][i] = 0;
			}
		}

		MainGraph[0][num_of_vertexes] = MainGraph[num_of_vertexes][0] = Vertexes[num_of_vertexes].Name = num_of_vertexes;

	}


	// DIJKSTRA'S ALGORITHM FUNCTIONS	
	//
	// Function complement to DijkstrasAl function that checks whether the element has already in the closed set
	bool check_path(int Vert)
	{
		int i = 1;

		while (Path[i]) {
			if (Path[i] == Vert)
				return 0;
			i++;
		}

		return 1;
	}

	// Recursive implementation of Dijkstra's algorithm - function
	void DijkstrasAl(int source)
	{
		Path[++path_counter] = source;						// When start - add a source vertex to the path and increase a path counter

		if (Vertexes[source].is_destination)				// If a source vertex is a destination vertex, THAN...
		{
			ShortestPath = Path;							// We have a one of the paths (perhabs the shortest, but need to check more)
			Path[path_counter--] = 0;						// Now we have stored the shortest path and need to change a path to a step behind
		}
		else {												// ELSE... 

			for (int i = 1; i <= num_of_vertexes; i++) {	// ...need to check all the vertexes in the region

				if ((MainGraph[source][i] != 0) && check_path(i)) {												// If the particular edge exists AND vertex i isn't in the closed set, than..

					if (Vertexes[source].Cost + MainGraph[source][i] < Vertexes[i].Cost) {							// If current cost for path to the checking vertex is greater than checking path than...
						Vertexes[i].Cost = Vertexes[source].Cost + MainGraph[source][i];							//...re-write a new cost to the checking vertex according to this path
						DijkstrasAl(i);																				// take this vertex and check its region
					}

				}
			}

			Path[path_counter--] = 0;																				// No more cheaper paths in the region - need to go to the root (source vertex of this path)

		}
	}


	// MST ALGORITHMS FUNCTIONS
	//
	// Returns zero graph with two vertexes
	Graph ReturnNullGraph()
	{
		Graph temp(0.0, 0, 2);

		return temp;

	}

	// Bubble sorting procedure for vector of graphs
	void bubble_sort(vector<Graph>& SimpleTrees)
	{
		Graph temp(0.0, 0, 2);

		for (int i = 0; i < SimpleTrees.size() - 1; i++)
			for (int j = SimpleTrees.size() - 1; j >= i + 1; j--)
				if (SimpleTrees[i].GetEdgeValue(1, 2) > SimpleTrees[j].GetEdgeValue(1, 2)) {
					temp = SimpleTrees[i];
					SimpleTrees[i] = SimpleTrees[j];
					SimpleTrees[j] = temp;
				}
	}

	// Function that indicates whether a loop in MST when try to add a new edge. Loop will be if we add X-Y edge to the MST and (1) there is X and Y edges in current tree and (2) there is path between X and Y
	bool check_for_loops(vector<Graph>& SimpleTrees, const int k, Graph& MST_Graph)
	{

		int is_First_Node = 0, is_Second_Node = 0, j = 0;										// Variables for testing whether X and Y already in a closed set

		while (((!is_First_Node) || (!is_Second_Node)) && (j < k)) {							// Check whether X and Y are in the closed set

			if ((SimpleTrees[k].GetVertexName(1) == SimpleTrees[j].GetVertexName(1))
				||
				(SimpleTrees[k].GetVertexName(1) == SimpleTrees[j].GetVertexName(2)))
				is_First_Node = 1;

			if ((SimpleTrees[k].GetVertexName(2) == SimpleTrees[j].GetVertexName(1))
				||
				(SimpleTrees[k].GetVertexName(2) == SimpleTrees[j].GetVertexName(2)))
				is_Second_Node = 1;

			j++;

		}

		if (is_First_Node && is_Second_Node) {													// If X and Y are in the closed set - try to find a path between X and Y using Dijkstra's algorithm

			if (MST_Graph.GetMinPathCost(SimpleTrees[k].GetVertexName(1), SimpleTrees[k].GetVertexName(2)))
				return 1;
			else
				return 0;
		}

		else
			return 0;

	}

	// Analyzing connections in initial graph. If one or more vertexes are disconnected - can't construct MST
	bool is_graph_disconnected()
	{

		for (int i = 2; i <= num_of_vertexes; i++)
		{
			if (!this->GetMinPathCost(1, i))
				return 1;
		}

		return 0;

	}



public:

	// CONSTRUCTORS OF CLASS GRAPH
	//
	// Constructor by default
	Graph()
	{
		create_zero_graph();
	}

	//Constructor with initial conditions from user: DENSITY - probability [0; 1] that particular vertex will have edges to all other vertexes. MAX_COST - maximum cost for path from one vertex to any other
	Graph(const double Density, const int MAX_COST)
	{
		srand(time(NULL));

		num_of_vertexes = MIN_NUM_OF_VERTEX + (rand() % (MAX_NUM_OF_VERTEX - MIN_NUM_OF_VERTEX + 1));							// Generate how many vertxes we will have (int from [MIN; MAX])

		MainGraph.resize(num_of_vertexes + 1, vector<int>(num_of_vertexes + 1));												// Creates graph in NxN matrix representation

		for (int i = 1; i < num_of_vertexes; i++) {
			MainGraph[0][i] = MainGraph[i][0] = Vertexes[i].Name = i;															// Generate the names of each vertex (from 1 to num_of_vertexes)

			for (int j = i + 1; j <= num_of_vertexes; j++) {
				if (((static_cast<double>(static_cast<int>(rand() % 100)) + 1) / 100) <= Density) {								// If there is no self-loop AND probability of edge <= Density THAN...
					MainGraph[i][j] = MainGraph[j][i] = rand() % MAX_COST + 1;													// Create edge and generate cost [1; MAX_COST from User]
					num_of_edges++;
				}
				else
					MainGraph[i][j] = MainGraph[j][i] = 0;																		// ELSE no edge - means cost is 0
			}
		}

		MainGraph[0][num_of_vertexes] = MainGraph[num_of_vertexes][0] = Vertexes[num_of_vertexes].Name = num_of_vertexes;

	}

	//Constructor with initial conditions from user: DENSITY - probability [0; 1] that particular vertex will have edges to all other vertexes. MAX_COST - maximum cost for path from one vertex to any other. EXACT_VERTEXES - how many vertexes will be in a graph 
	Graph(const double Density, const int MAX_COST, const int exact_vertexes)
	{
		if ((exact_vertexes > MAX_NUM_OF_VERTEX) || (exact_vertexes < MIN_NUM_OF_VERTEX)) {
			cout << "Can't create a graph. Number of vertexes should be within [" << MIN_NUM_OF_VERTEX << " ; " << MAX_NUM_OF_VERTEX << "]" << endl;
			create_zero_graph();
			cout << "Zero cost graph with 2 vertexes has been created!" << endl;

		}
		else {

			num_of_vertexes = exact_vertexes;																					// How many vertexes will be in a Graph, explicitly determined by User

			MainGraph.resize(num_of_vertexes + 1, vector<int>(num_of_vertexes + 1));											// Creates graph in NxN matrix representation

			for (int i = 1; i < num_of_vertexes; i++) {
				MainGraph[0][i] = MainGraph[i][0] = Vertexes[i].Name = i;														// Generate the names of each vertex (from 1 to num_of_vertexes)

				for (int j = i + 1; j <= num_of_vertexes; j++) {
					if (((static_cast<double>(static_cast<int>(rand() % 100)) + 1) / 100) <= Density) {							// If probability of edge existing <= Density THAN...
						MainGraph[i][j] = MainGraph[j][i] = rand() % MAX_COST + 1;												// Creates edge and generates cost for that [1; MAX_COST from User]
						num_of_edges++;																							// Increase number of edges
					}
					else
						MainGraph[i][j] = MainGraph[j][i] = 0;																	// ELSE no edge - means cost is 0
				}
			}

			MainGraph[0][num_of_vertexes] = MainGraph[num_of_vertexes][0] = Vertexes[num_of_vertexes].Name = num_of_vertexes;

		}
	}

	// Constructor for reading Graph from file "File_name"
	Graph(string File_name)
	{
		// Trying to open a file.
		ifstream Input_File(File_name, ios_base::in);

		if (Input_File.fail()) {
			cout << "Can't open file " << File_name << " . Please check the path to the file." << endl;
			create_zero_graph();
			cout << "Zero cost graph with 2 vertexes has been created!" << endl;

		}
		else {

			istream_iterator <int> start(Input_File), end;
			vector <int> Data(start, end);

			num_of_vertexes = Data[0];

			MainGraph.resize(num_of_vertexes + 1, vector<int>(num_of_vertexes + 1));											// Creates graph in NxN matrix representation

			for (int i = 1; i <= num_of_vertexes; i++)																			// Give names to all vertixes in Graph
				MainGraph[0][i] = MainGraph[i][0] = Vertexes[i].Name = i;

			for (int i = 1; i < Data.size(); i += 3)
				if (!this->GetEdgeValue(Data[i] + 1, Data[i + 1] + 1))
					this->SetEdgeValue(Data[i] + 1, Data[i + 1] + 1, Data[i + 2]);												// "+1" means that this class starting calcs in Graph from vertex 1, not from 0

		}
	}


	//OVERLOADING OPERATORS FOR CLASS GRAPH
	//
	// Overloading << operator - print out a whole graph on screen. IMPORTANT: can print-out a graph in incorrected way (depends on system). For full correctness output please use method "CopyToFile" and see the graph to the txt file
	friend ostream& operator<<(ostream& out, Graph& Gr1)
	{
		for (int i = 0; i <= Gr1.GetVertexes(); i++) {
			for (int j = 0; j <= Gr1.GetVertexes(); j++) {

				if (Gr1.GetEdgeValue(i, j))
					cout << Gr1.GetEdgeValue(i, j) << "\t";
				else
					cout << "-\t";
			}

			cout << endl;
		}

		cout << "Num. of vertexes : " << Gr1.GetVertexes() << endl;
		cout << "Num. of edges: " << Gr1.GetEdges() << endl;
		cout << "Total graph cost: " << Gr1.GetTotalGraphCost() << endl;

		return out;
	}

	// Overloading = operator
	Graph& operator=(Graph& Gr1)
	{
		num_of_vertexes = Gr1.GetVertexes();
		num_of_edges = Gr1.GetEdges();
		Vertex* Vert_temp = Gr1.GetVertexesContent();

		// Resize initial matrix
		MainGraph.clear();
		MainGraph.resize(num_of_vertexes + 1, vector<int>(num_of_vertexes + 1));

		//Copying all edges
		for (int i = 0; i <= num_of_vertexes; i++)
			for (int j = 0; j <= num_of_vertexes; j++)
				MainGraph[i][j] = Gr1.GetEdgeValue(i, j);

		//Copying vertex info
		for (int i = 0; i <= MAX_NUM_OF_VERTEX; i++) {
			Vertexes[i].Name = Vert_temp[i].Name;
			Vertexes[i].Cost = Vert_temp[i].Cost;
			Vertexes[i].is_destination = Vert_temp[i].is_destination;
		}

		return *this;
	}


	// PUBLIC GENERAL METHODS OF CLASS GRAPH
	//
	// Copy Graph to a file for analysis and checking results
	void CopyToFile(string FileName)
	{
		ofstream Output_file(FileName, ios_base::out);

		for (int i = 0; i <= num_of_vertexes; i++) {
			for (int j = 0; j <= num_of_vertexes; j++)
				Output_file << MainGraph[i][j] << "\t";

			Output_file << endl;
		}

		Output_file << "Num. of vertexes: " << num_of_vertexes << endl;
		Output_file << "Num. of edges: " << num_of_edges << endl;
		Output_file << "Total Graph cost: " << this->GetTotalGraphCost() << endl;

		cout << "Graph has been copied to a file: " << FileName << endl;
	}

	// Returns a number of edges in a graph
	int GetEdges() { return num_of_edges; }

	// Returns a num of vertexes in a graph
	int GetVertexes() { return num_of_vertexes; }

	// Returns cost value of X to Y edge
	int GetEdgeValue(const int x, const int y) { return MainGraph[x][y]; }

	// Returns the name of X
	int GetVertexName(const int x) { return MainGraph[0][x]; }

	// Returns a pointer to Vertex object 
	Vertex* GetVertexesContent() { return Vertexes; }

	// Returns bool whether edge between X and Y exists
	bool GetAdjacent(const int x, const int y)
	{
		if (MainGraph[x][y])
			return 1;
		else
			return 0;
	}

	// Returns sum of costs of all edges in a graph
	int GetTotalGraphCost()
	{
		int sum = 0;
		for (int i = 1; i < num_of_vertexes; i++)
			for (int j = i + 1; j <= num_of_vertexes; j++)
				sum += MainGraph[i][j];

		return sum;
	}

	// Prints out list of all neighborhood vertexes to X
	void Neighbors(const int x)
	{
		cout << "All neighborhoods to vertex " << x << ":" << endl;
		for (int i = 1; i <= num_of_vertexes; i++) {
			if (MainGraph[x][i])
				cout << "Vertex " << i << " , cost: " << MainGraph[x][i] << endl;
		}
	}

	// Adds X to Y edge of value COST if this edge doesn't exist
	void AddEdge(const int x, const int y, const int cost)
	{
		if (!MainGraph[x][y]) {
			MainGraph[x][y] = MainGraph[y][x] = cost;
			num_of_edges++;
		}
		else
			cout << "The edge between " << x << " and " << y << " has already exists. Can't add it." << endl;
	}

	// Deletes edge between X to Y if this edge exists
	void Delete(const int x, const int y)
	{
		if (MainGraph[x][y]) {
			MainGraph[x][y] = MainGraph[y][x] = 0;
			num_of_edges--;
		}
		else
			cout << "The edge between " << x << " and " << y << "doesn't exists. Can't delete the edge." << endl;
	}

	// Deteles all edges for particular vertex (V) in a graph
	void DeleteEdges(const int& V)
	{
		for ( int i = 1; i <= num_of_vertexes; i++ )
			if (MainGraph[V][i] != 0) {
				MainGraph[V][i] = MainGraph[i][V] = 0;
				num_of_edges--;
			}
	}


	// Sets X to Y edge value 
	void SetEdgeValue(const int x, const int y, const int Val)
	{
		if ((x <= num_of_vertexes) && (y <= num_of_vertexes)) {
			MainGraph[x][y] = MainGraph[y][x] = Val;
			num_of_edges++;
		}
		else
			cout << "Can't set edge between " << x << " and " << y << " because one of the vertexes doesn't exist." << endl;
	}

	// Sets the name (int) for particular node
	void SetName(int node, int Name)
	{
		MainGraph[node][0] = MainGraph[0][node] = Name;
	}

	// Deletes all edges in a graph
	void DeleteAllEdges()
	{
		for (int i = 1; i < num_of_vertexes; i++)
			for (int j = i + 1; j <= num_of_vertexes; j++) {
				if (MainGraph[i][j]) {
					MainGraph[i][j] = MainGraph[j][i] = 0;
					num_of_edges--;
				}
			}
	}

	// Make a resize zero (empty) graph to N vertexes
	void Resize(int N)
	{
		if ((N < MIN_NUM_OF_VERTEX) || (N > MAX_NUM_OF_VERTEX))
			cout << "Can't make a resize, because number of vertexes should be within [" << MIN_NUM_OF_VERTEX << " ; " << MAX_NUM_OF_VERTEX << "]" << endl;
		else if (num_of_edges)
			cout << "Can't make a resize, because graph is not a zero graph. First make a DeleteAllEdges method." << endl;
		else {
			num_of_vertexes = N;

			MainGraph.clear();
			MainGraph.resize(num_of_vertexes + 1, vector<int>(num_of_vertexes + 1));

			for (int i = 1; i < num_of_vertexes; i++) {
				MainGraph[0][i] = MainGraph[i][0] = Vertexes[i].Name = i;

				for (int j = i + 1; j <= num_of_vertexes; j++)
					MainGraph[i][j] = MainGraph[j][i] = 0;
			}

			MainGraph[0][num_of_vertexes] = MainGraph[num_of_vertexes][0] = Vertexes[num_of_vertexes].Name = num_of_vertexes;

		}

	}

	// Adds new vertex into graph to the end 
	void AddVertex()
	{

		vector <vector <int>> temp_graph = MainGraph;

		num_of_vertexes++;

		MainGraph.clear();
		MainGraph.resize(num_of_vertexes + 1, vector<int>(num_of_vertexes + 1));

		for (int i = 1; i < num_of_vertexes; i++) {
			MainGraph[0][i] = MainGraph[i][0] = Vertexes[i].Name = i;														// Generate the names of each vertex (from 1 to num_of_vertexes)

			for (int j = i + 1; j <= num_of_vertexes; j++) {
				if ((i != num_of_vertexes) && (j != num_of_vertexes))
					MainGraph[i][j] = MainGraph[j][i] = temp_graph[i][j];
			}
		}

		MainGraph[0][num_of_vertexes] = MainGraph[num_of_vertexes][0] = Vertexes[num_of_vertexes].Name = num_of_vertexes;

	}

	
	// METHODS FOR DIJKSTRA'S SHOREST PATH ALGORITHM
	// 
	// Finds and prints out the shortest path in the graph between SOURCE and DESTINATION vertexes using Dijkstra's algorithm
	string GetShortestPath(const int source, const int destination)
	{
		// Preparations for calcs
		string return_string;
		Vertexes[source].Cost = 0;
		Vertexes[destination].is_destination = 1;
		Path.resize(num_of_vertexes + 1, 0);						// Change vector size for testing set
		ShortestPath.resize(num_of_vertexes + 1, 0);				// Change vector size for final closed set

		DijkstrasAl(source);										// Function DijkstraAl provides calculations of shortest path and minimum cost using Dijkstra's Algorithm

		if (Vertexes[destination].Cost != Vertexes[destination].GetInitCost()) {

			return_string = to_string(ShortestPath[1]);

			int i = 2;
			while (ShortestPath[i])
				return_string = return_string + " - " + to_string(ShortestPath[i++]);
			cout << endl;
		}
		else
			return_string = "No such path.";


		// Clean up		
		for (int i = 1; i <= num_of_vertexes; i++)
			Vertexes[i].Cost = Vertexes[i].GetInitCost();

		Vertexes[destination].is_destination = 0;

		return return_string;

	}

	// Finds and returns the smallest cost of path between SOURCE and DESTINATION vertexes in the graph using Dijkstra's algorithm
	int GetMinPathCost(const int source, const int destination)
	{
		// Preparations for calcs
		Vertexes[source].Cost = 0;
		Vertexes[destination].is_destination = 1;
		Path.resize(num_of_vertexes + 1, 0);						// Vector for testing set
		ShortestPath.resize(num_of_vertexes + 1, 0);				// Vector for final closed set

		DijkstrasAl(source);										// Function DijkstraAl provides calculations of shortest path and minimum cost using Dijkstra's Algorithm

		if (Vertexes[destination].Cost != Vertexes[destination].GetInitCost()) {

			int Real_cost = Vertexes[destination].Cost;
			// Clean up before returning
			for (int i = 1; i <= num_of_vertexes; i++)
				Vertexes[i].Cost = Vertexes[i].GetInitCost();

			Vertexes[destination].is_destination = 0;

			return Real_cost;
		}
		else {

			// Clean up before returning
			for (int i = 1; i <= num_of_vertexes; i++)
				Vertexes[i].Cost = Vertexes[i].GetInitCost();

			Vertexes[destination].is_destination = 0;

			return 0;
		}

	}

	// Special Dijkstra's algorithm for HEX game finding a winner
	int SpecialHexAlgorithm(vector <int> sources, vector <int> destinations)
	{
		// Preparations for calcs
		Path.resize(num_of_vertexes + 1, 0);						// Vector for testing set
		ShortestPath.resize(num_of_vertexes + 1, 0);				// Vector for final closed set

		// Make all destinations active
		for (int i = 0; i < destinations.size(); i++) 
			Vertexes[destinations[i]].is_destination = 1;

		// Start to calculate shortest path 
		for (int i = 0; i < sources.size(); i++) {

			Vertexes[sources[i]].Cost = 0;

			DijkstrasAl(sources[i]);				// Function DijkstraAl provides calculations of shortest path and minimum cost using Dijkstra's Algorithm

			// Check whether any shortest path from source to any destinations
			for (int i = 0; i < destinations.size(); i++) {
				if (Vertexes[destinations[i]].Cost != Vertexes[destinations[i]].GetInitCost()) {

					int Real_cost = Vertexes[destinations[i]].Cost;

					// Clean up before returning
					for (int i = 1; i <= num_of_vertexes; i++) 
						Vertexes[i].Cost = Vertexes[i].GetInitCost();
						
					for (int i = 0; i < destinations.size(); i++) {
						Vertexes[destinations[i]].is_destination = 0;
						
					}

					return Real_cost;
				}
			}

		}

		// If we are here - it means there is no shortest pathes among sources and destsinations. Clean up before returning
		for (int i = 1; i <= num_of_vertexes; i++) {
			Vertexes[i].Cost = Vertexes[i].GetInitCost();
			Vertexes[i].is_destination = 0;
		}

		return 0;
			
	}


	// METHODS FOR KRUSKAL'S MST ALGORITHM
	// 
	Graph& CreateMST()
	{

		if (is_graph_disconnected()) {
			cout << "Initial graph is disconnected! Can't create MST form it." << endl;
			//return this->ReturnNullGraph();														// Just returns a zero graph with two vertexes
		}
		else {

			int count = 0;
			num_of_simple_trees = num_of_edges;
			SimpleTrees.resize(num_of_simple_trees);

			Graph Temp_Graph;																	// Graph for temporary storing current tree (edge)

			for (int i = 1; i < num_of_vertexes; i++)											// Creates all possible simple trees from a graph
				for (int j = i + 1; j <= num_of_vertexes; j++)

					if (MainGraph[i][j]) {

						Temp_Graph.SetName(1, MainGraph[i][0]);
						Temp_Graph.SetName(2, MainGraph[j][0]);
						Temp_Graph.SetEdgeValue(1, 2, MainGraph[i][j]);

						SimpleTrees[count++] = Temp_Graph;										// All simple trees are storing in SimpleTrees vector

					}

			bubble_sort(SimpleTrees);															// Make a sorting procedure for all trees by its edge value (cost)


			static Graph MST(0.0, 0, num_of_vertexes);											// Create non-distructable variable for returning MST from class Graph
			MST.DeleteAllEdges();																// Make sure we have empty zero graph...
			MST.Resize(num_of_vertexes);														//...with exact num. of vertexes

			// Trying to construct MST from sorted simple trees vector
			for (int i = 0; i < num_of_simple_trees; i++)
				if ((i == 0) || (i == 1)) {														// Granted cases that there are no loops
					MST.SetEdgeValue(
						SimpleTrees[i].GetVertexName(1),
						SimpleTrees[i].GetVertexName(2),
						SimpleTrees[i].GetEdgeValue(1, 2));

				}
				else {
					if (!check_for_loops(SimpleTrees, i, MST)) {

						MST.SetEdgeValue(
							SimpleTrees[i].GetVertexName(1),
							SimpleTrees[i].GetVertexName(2),
							SimpleTrees[i].GetEdgeValue(1, 2));

					}
				}

			return MST;																			// Returns constructing MST graph

		}

	}


	// DESTRUCTOR for cleaning heap
	~Graph()
	{
		delete[] Vertexes;
	}

};
