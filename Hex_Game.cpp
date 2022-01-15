#include <iostream>
#include <vector>
#include <string>
#include "Graph.h"

using namespace std;

//CLASS HEXVERTEX. Needs for storing information about all nodes on the hexplot
class HexVertex {

private:

	int NodeName,		// The name (number) of the node from 1 to size^2. NodeName = (x-1)*size + y
		ElemNum_in_graph, // Stores the number of the node in main graph for quick searching
		x,				// X-coorinate of the node on the hexplot
		y,				// Y-coorinate of the node on the hexplot
		PlotSize,		// Size of hexplot
		Owner;			// Who owns the node (what player: "B" - 1 or "W" - 2 or nobody - 0)


public:

	HexVertex()
	{
		NodeName = 0;
		Owner = '-';
		x = y = 0;
		PlotSize = 0;
		ElemNum_in_graph = 0;
	}

	int GetName() { return NodeName; }

	int GetX() { return x; }

	int GetY() { return y; }

	int GetOwner() { return Owner; }

	// Returns the number of particular vertex in graph. Need for DIjkstra's algorithm to find a path and to know who is a winner
	int GetNumOfElementInGraph() { return ElemNum_in_graph; }

	// Idea is to determine the name as a number of particular node on hexplot, where (1,1) is number 1 and (size, size) is number size*size
	void SetXYandName(const int& New_X, const int& New_Y, const int& Size)
	{
		x = New_X;
		y = New_Y;
		NodeName = (x - 1) * Size + y;
		PlotSize = Size;
	}

	void SetOwner(int Ow) { Owner = Ow; }

	// Sets the number of particular vertex in graph. Need for DIjkstra's algorithm to find a path and to know who is a winner
	void SetNumOfElementInGraph(const int& elem_num) { ElemNum_in_graph = elem_num; }

};


class HexPlot : public Graph
{

private:

	// Private variables of class HexPlot
	int max_num_of_nodes,			// Maximum theoretical amount of activated nodes in the game
		active_num_of_nodes,		// Current activated (has an owner) nodes on the hexplot		
		PlotSize;					// The size of a hexplot


	vector <HexVertex> Nodes;		// Contains info about all active nodes on hexplot
	int Max_element;				// Sotres the value of a name (=number) of maximum element in the graph

	// Function checks whether X-Y node is empty or not
	int is_empty_place(const int& Check_X, const int& Check_Y)
	{
		if (((Check_X <= PlotSize) && (Check_Y <= PlotSize)) && ((Check_X >= 1) && (Check_Y >= 1)))

			if (Nodes[(Check_X - 1) * PlotSize + Check_Y].GetName() != 0)
				return 0;				// Node is NOT empty
			else
				return 1;				// Node is empty
	}


public:

	// Constructor by defalt - creates hexplot for a game size 2x2
	HexPlot()
	{
		PlotSize = 2;
		max_num_of_nodes = PlotSize * PlotSize;
		active_num_of_nodes = 0;
		Nodes.clear();
		Nodes.resize(401);					// Max size of 20x20 hexplot	
		Max_element = 0;

		// Zero element (Node[0]) as a null analogue					

		// Min empty Graph 2x2 has been created by base class Graph constructor
	}

	// Constructor that creates hexplot for a game (size X size)
	HexPlot(const int& size)
	{
		PlotSize = size;
		max_num_of_nodes = PlotSize * PlotSize;
		active_num_of_nodes = 0;
		Nodes.clear();
		Nodes.resize(401);					// Max size of 20x20 hexplot
		Max_element = 0;

		// Zero element (Node[0]) as a null analogue

		// Min empty Graph 2x2 has been created by base class Graph constructor

	}

	//GENERAL METHODS OF CLASS HEXPLOT
	//
	// Returns number of current activated nodes on the hexplot of all players
	int GetNumOfActiveNodes() { return active_num_of_nodes; }

	// Returns number of maximum nodes on hexplot
	int GetMaxNumOfNodes() { return max_num_of_nodes; }

	// Returns current maximum element number on the hexplot
	int GetCurrentMaxEl() { return Max_element; }

	// Returns the size of hexplot
	int GetPlotSize() { return PlotSize; }

	// Returns a pointer to a current Node using X-Y coordinates from hexplot
	HexVertex* GetNodeInfoByXY(const int& x, const int& y) { return &Nodes[(x - 1) * PlotSize + y]; }

	// Set the initial position for nodes vector
	void DeleteAllNodes()
	{
		Nodes.clear();
		Nodes.resize(401);
	}

	// Change a hexplot size. IMPORTANT! Method erase all info about current game reflected to hexplot
	void SetPlotSize(const int& New_Size)
	{
		PlotSize = New_Size;
		max_num_of_nodes = PlotSize * PlotSize;
		active_num_of_nodes = 0;
		Nodes.clear();
		Nodes.resize(401);					// Max size of 20x20 hexplot

		// Total clearing of a graph
		this->DeleteAllEdges();
		this->Resize(2);
	}

	// Method that includes a new node into game after any player's move
	int NewMove(const int& New_X, const int& New_Y, const int& Pl_color)
	{

		// Checking whether X-Y is a legal move
		if (((New_X <= PlotSize) && (New_Y <= PlotSize)) && ((New_X >= 1) && (New_Y >= 1)) && is_empty_place(New_X, New_Y)) {

			// Set values for a new node in Nodes vector
			Nodes[(New_X - 1) * PlotSize + New_Y].SetXYandName(New_X, New_Y, PlotSize);
			Nodes[(New_X - 1) * PlotSize + New_Y].SetOwner(Pl_color);
			Nodes[(New_X - 1) * PlotSize + New_Y].SetNumOfElementInGraph(++active_num_of_nodes);

			// Check for a new Max_element in Nodes vector
			if (Max_element < ((New_X - 1) * PlotSize + New_Y))
				Max_element = (New_X - 1) * PlotSize + New_Y;

			// Create this node in the graph
			if (active_num_of_nodes > 2)
				this->AddVertex();					// We have 2 nodes at initialization. If one more node is coming - need to increase number of nodes in the graph

			// Check all possible connections of a new node and create new edges in the graph
			if ((New_X != 1) && Nodes[((New_X - 1) - 1) * PlotSize + New_Y].GetNumOfElementInGraph() && (Nodes[((New_X - 1) - 1) * PlotSize + New_Y].GetOwner() == Pl_color))		// Case for neighborhood node with coordinates (x-1, y) AND neighborhood is exists AND player's color == neighborhood's color
				this->AddEdge(Nodes[((New_X - 1) - 1) * PlotSize + New_Y].GetNumOfElementInGraph(), active_num_of_nodes, Pl_color);			// Make an edge in a graph between these nodes for cost "Pl_color"			

			if ((New_Y != 1) && Nodes[(New_X - 1) * PlotSize + (New_Y - 1)].GetNumOfElementInGraph() && (Nodes[(New_X - 1) * PlotSize + (New_Y - 1)].GetOwner() == Pl_color)) 		// Case for neighborhood node with coordinates (x, y-1) AND neighborhood is exists AND player's color == neighborhood's color
				this->AddEdge(Nodes[(New_X - 1) * PlotSize + (New_Y - 1)].GetNumOfElementInGraph(), active_num_of_nodes, Pl_color);			// Make an edge in a graph between these nodes for cost "Pl_color"


			if ((New_X != PlotSize) && Nodes[((New_X - 1) + 1) * PlotSize + New_Y].GetNumOfElementInGraph() && (Nodes[((New_X - 1) + 1) * PlotSize + New_Y].GetOwner() == Pl_color)) 		// Case for neighborhood node with coordinates (x+1, y) AND neighborhood is exists AND player's color == neighborhood's color
				this->AddEdge(Nodes[((New_X - 1) + 1) * PlotSize + New_Y].GetNumOfElementInGraph(), active_num_of_nodes, Pl_color);			// Make an edge in a graph between these nodes for cost "Pl_color"


			if ((New_Y != PlotSize) && Nodes[(New_X - 1) * PlotSize + (New_Y + 1)].GetNumOfElementInGraph() && (Nodes[(New_X - 1) * PlotSize + (New_Y + 1)].GetOwner() == Pl_color)) 		// Case for neighborhood node with coordinates (x, y+1) AND neighborhood is exists AND player's color == neighborhood's color
				this->AddEdge(Nodes[(New_X - 1) * PlotSize + (New_Y + 1)].GetNumOfElementInGraph(), active_num_of_nodes, Pl_color);			// Make an edge in a graph between these nodes for cost "Pl_color"


			if ((New_X != 1) && (New_Y != PlotSize) && Nodes[((New_X - 1) - 1) * PlotSize + (New_Y + 1)].GetNumOfElementInGraph() && (Nodes[((New_X - 1) - 1) * PlotSize + (New_Y + 1)].GetOwner() == Pl_color))			// Case for neighborhood node with coordinates (x-1, y+1) AND neighborhood is exists AND player's color == neighborhood's color
				this->AddEdge(Nodes[((New_X - 1) - 1) * PlotSize + (New_Y + 1)].GetNumOfElementInGraph(), active_num_of_nodes, Pl_color);			// Make an edge in a graph between these nodes for cost "Pl_color"


			if ((New_X != PlotSize) && (New_Y != 1) && Nodes[((New_X - 1) + 1) * PlotSize + (New_Y - 1)].GetNumOfElementInGraph() && (Nodes[((New_X - 1) + 1) * PlotSize + (New_Y - 1)].GetOwner() == Pl_color))			// Case for neighborhood node with coordinates (x+1, y-1) AND neighborhood is exists AND player's color == neighborhood's color
				this->AddEdge(Nodes[((New_X - 1) + 1) * PlotSize + (New_Y - 1)].GetNumOfElementInGraph(), active_num_of_nodes, Pl_color);		// Make an edge in a graph between these nodes for cost "Pl_color"

			return 1;
		}
		else {
			if (is_empty_place(New_X, New_Y)) {
				cout << "Can't add a new node because X or Y is out of the hexplot size. Please make another move" << endl;
				return -1;
			}
			else {
				cout << "Node with coordinates (" << New_X << " ; " << New_Y << ") has already exists. Please make another move" << endl;
				return 0;
			}
		}

	}

	// Method that returns the status of winnig taking into account the last move
	int isWin(const int& New_X, const int& New_Y, const int& Pl_color)
	{
		int Win = 0, i = 1, j = 1;

		if (Pl_color == 1) {			// Check for winning for N-S Player

			while ((i <= PlotSize) && !Win) {  //Check all N-nodes

				if (Nodes[i].GetName()) {		// If N-node exists

					while ((j <= PlotSize) && !Win) {		// Check all S-nodes 

						if (Nodes[(PlotSize - 1) * PlotSize + j].GetName())		// If any S-node is exists too THAN...
							Win = this->GetMinPathCost(Nodes[i].GetNumOfElementInGraph(), Nodes[(PlotSize - 1) * PlotSize + j].GetNumOfElementInGraph());		//...try to find a path between N and S nodes. If so - Player1 is a winner

						j++;
					}

				}

				i++;
				j = 1;

			}

		}
		else {							// Check for winning for W-E Player

			while ((i <= PlotSize) && !Win) {		//Check all W-nodes

				if (Nodes[(i - 1) * PlotSize + 1].GetName()) {		// If W-node is exists

					while ((j <= PlotSize) && !Win) {		// Check all E-nodes 

						if (Nodes[(j - 1) * PlotSize + PlotSize].GetName())		// If any E-node is exists TNAN...
							Win = this->GetMinPathCost(Nodes[(i - 1) * PlotSize + 1].GetNumOfElementInGraph(), Nodes[(j - 1) * PlotSize + PlotSize].GetNumOfElementInGraph());		//...try to find a path between W and E nodes. If so - Player2 is a winner

						j++;
					}

				}

				i++;
				j = 1;

			}

		}

		return Win;

	}

	// Generate random filled hexplot according to current hexplot position
	void MakeRandomDistr()
	{
		int what_color,															// Who makes a move (1 - Black, 2 - White)
			count_W = (max_num_of_nodes - active_num_of_nodes) / 2,				// How many W moves remain to the end of game
			count_B = (max_num_of_nodes - active_num_of_nodes) - count_W;		// How many B moves remain to the end of game

		for (int i = 1; i <= PlotSize; i++)
			for (int j = 1; j <= PlotSize; j++)
				if (!Nodes[(i - 1) * PlotSize + j].GetName()) {

					if (count_W != 0 and count_B != 0) {

						what_color = rand() % 2 + 1;			// Distribution function. Can take two values: 1 or 2. I'm thinking about its improvement...
						
						this->NewMove(i, j, what_color);		// Make a random move (shuffling)

						if (what_color == 2)
							count_W--;
						else
							count_B--;

					}
					else if (!count_W) {
						this->NewMove(i, j, 1);
						count_B--;
					}
					else if (!count_B) {
						this->NewMove(i, j, 2);
						count_W--;
					}

				}
	}

	// Check who wins after current move
	int WhoWins()
	{
		int Win = 0, i = 1, j = 1, c = 0,
			num_of_sources = 0,					// How many North or West nodes exist
			num_of_destinations = 0;			// How many South or East nodes exist
		vector <int> Sources_N(122), Sources_W(122), Destinations_S(122), Destinations_E(122);						// Vector with all available destinations for further special Dijkstra's algorithm

		// Check for winning for N-S Player (Player1 or Black-player)
		// Create vector of all Source nodes (all North nodes owned by Player1)
		for (int i = 1; i <= PlotSize; i++)
			if (Nodes[i].GetName() && Nodes[i].GetOwner() == 1) {				// If particular North-node for Player1 exists
				Sources_N[num_of_sources] = Nodes[i].GetNumOfElementInGraph();					// Store such node (its number in graph) in "Sources" vector
				num_of_sources++;
			}
		
		if (num_of_sources) {												// If we have one or more sources, THAN let's create vector of all destination nodes (South-nodes) that owned by Player1

			for (int i = (PlotSize - 1) * PlotSize + 1; i <= PlotSize * PlotSize; i++)
				if (Nodes[i].GetName() && Nodes[i].GetOwner() == 1) {				// If particular South-node for Player1 exists
					Destinations_S[num_of_destinations] = Nodes[i].GetNumOfElementInGraph();					// Store such node (its number in graph) in "Destinations" vector
					num_of_destinations++;
				}
			
			if (num_of_destinations)									// If we also has one or more destinations, THAN... 
				Win = this->SpecialHexAlgorithm(Sources_N, Destinations_S);		// ...check the path among all Sources and Destinations
		}
			

		// if Player1 is a winner THAN stop function and return Win (value is 1) to the parent function
		if (Win)
			return Win = 1;


		// If NOT - check for winning for W-E Player (Player2 or White-player)
		num_of_sources = num_of_destinations = 0;

		// Create vector of all Source nodes (all West nodes owned by Player2)
		for (int i = 1; i <= PlotSize * PlotSize; i += PlotSize)
			
			if (Nodes[i].GetName() && Nodes[i].GetOwner() == 2) {				// If particular West-node for Player2 exists
				Sources_W[num_of_sources] = Nodes[i].GetNumOfElementInGraph();					// Store such node (its number in graph) in "Sources" vector
				num_of_sources++;
			}
		
		if (num_of_sources) {												// If we have one or more sources, THAN lte's create vector of all destination nodes (West-nodes) that owned by Player2

			for (int i = PlotSize; i <= PlotSize * PlotSize; i += PlotSize)
				if (Nodes[i].GetName() && Nodes[i].GetOwner() == 2) {				// If particular East-node for Player2 exists
					Destinations_E[num_of_destinations] = Nodes[i].GetNumOfElementInGraph();					// Store such node (its number in graph) in "Destinations" vector
					num_of_destinations++;
				}
						
			if (num_of_destinations)									// If we also has one or more destinations, THAN... 
				Win = this->SpecialHexAlgorithm(Sources_W, Destinations_E);		// ...check the path among all Sources and DEstinations
		}
		

		// if Player2 is a winner THAN stop function and return Win (value is 2) to the parent function
		if (Win)
			return Win = 2;

		// If NOBODY wins THAN check who did the last move - it's a winner if no other moves on hexplot
		if (active_num_of_nodes == max_num_of_nodes)
			if (max_num_of_nodes % 2)
				return Win = 1;
			else
				return Win = 2;
		else
			return Win = 0;
	}


	// Cretaes all possible edges for particular node on hexplot's graph
	void CreateNewEdges(const int& Node_name)
	{				
		int
			New_X = Nodes[Node_name].GetX(),
			New_Y = Nodes[Node_name].GetY(),
			Pl_color = Nodes[Node_name].GetOwner();
								
		// Creation of new edges for particular node		
		if ((New_X != 1) && Nodes[((New_X - 1) - 1) * PlotSize + New_Y].GetNumOfElementInGraph() && !this->GetEdgeValue(Nodes[((New_X - 1) - 1) * PlotSize + New_Y].GetNumOfElementInGraph(), Nodes[Node_name].GetNumOfElementInGraph()) && (Nodes[((New_X - 1) - 1) * PlotSize + New_Y].GetOwner() == Pl_color)) 	// Case for neighborhood node with coordinates (x-1, y) AND neighborhood is exists AND edge doesn't exist yet AND player's color == neighborhood's color
			this->AddEdge(Nodes[((New_X - 1) - 1) * PlotSize + New_Y].GetNumOfElementInGraph(), Nodes[Node_name].GetNumOfElementInGraph(), Pl_color);			// Make an edge in a graph between these nodes for cost "Pl_color"			

		if ((New_Y != 1) && Nodes[(New_X - 1) * PlotSize + (New_Y - 1)].GetNumOfElementInGraph() && !this->GetEdgeValue(Nodes[(New_X - 1) * PlotSize + (New_Y - 1)].GetNumOfElementInGraph(), Nodes[Node_name].GetNumOfElementInGraph()) && (Nodes[(New_X - 1) * PlotSize + (New_Y - 1)].GetOwner() == Pl_color)) 	// Case for neighborhood node with coordinates (x, y-1) AND neighborhood is exists AND edge doesn't exist yet AND player's color == neighborhood's color
			this->AddEdge(Nodes[(New_X - 1) * PlotSize + (New_Y - 1)].GetNumOfElementInGraph(), Nodes[Node_name].GetNumOfElementInGraph(), Pl_color);			// Make an edge in a graph between these nodes for cost "Pl_color"

		if ((New_X != PlotSize) && Nodes[((New_X - 1) + 1) * PlotSize + New_Y].GetNumOfElementInGraph() && !this->GetEdgeValue(Nodes[((New_X - 1) + 1) * PlotSize + New_Y].GetNumOfElementInGraph(), Nodes[Node_name].GetNumOfElementInGraph()) && (Nodes[((New_X - 1) + 1) * PlotSize + New_Y].GetOwner() == Pl_color)) 		// Case for neighborhood node with coordinates (x+1, y) AND neighborhood is exists AND edge doesn't exist yet AND player's color == neighborhood's color
			this->AddEdge(Nodes[((New_X - 1) + 1) * PlotSize + New_Y].GetNumOfElementInGraph(), Nodes[Node_name].GetNumOfElementInGraph(), Pl_color);			// Make an edge in a graph between these nodes for cost "Pl_color"

		if ((New_Y != PlotSize) && Nodes[(New_X - 1) * PlotSize + (New_Y + 1)].GetNumOfElementInGraph() && !this->GetEdgeValue(Nodes[(New_X - 1) * PlotSize + (New_Y + 1)].GetNumOfElementInGraph(), Nodes[Node_name].GetNumOfElementInGraph()) && (Nodes[(New_X - 1) * PlotSize + (New_Y + 1)].GetOwner() == Pl_color)) 		// Case for neighborhood node with coordinates (x, y+1) AND neighborhood is exists AND edge doesn't exist yet AND player's color == neighborhood's color
			this->AddEdge(Nodes[(New_X - 1) * PlotSize + (New_Y + 1)].GetNumOfElementInGraph(), Nodes[Node_name].GetNumOfElementInGraph(), Pl_color);			// Make an edge in a graph between these nodes for cost "Pl_color"

		if ((New_X != 1) && (New_Y != PlotSize) && Nodes[((New_X - 1) - 1) * PlotSize + (New_Y + 1)].GetNumOfElementInGraph() && !this->GetEdgeValue(Nodes[((New_X - 1) - 1) * PlotSize + (New_Y + 1)].GetNumOfElementInGraph(), Nodes[Node_name].GetNumOfElementInGraph()) && (Nodes[((New_X - 1) - 1) * PlotSize + (New_Y + 1)].GetOwner() == Pl_color))			// Case for neighborhood node with coordinates (x-1, y+1) AND neighborhood is exists AND edge doesn't exist yet AND player's color == neighborhood's color
			this->AddEdge(Nodes[((New_X - 1) - 1) * PlotSize + (New_Y + 1)].GetNumOfElementInGraph(), Nodes[Node_name].GetNumOfElementInGraph(), Pl_color);			// Make an edge in a graph between these nodes for cost "Pl_color"

		if ((New_X != PlotSize) && (New_Y != 1) && Nodes[((New_X - 1) + 1) * PlotSize + (New_Y - 1)].GetNumOfElementInGraph() && !this->GetEdgeValue(Nodes[((New_X - 1) + 1) * PlotSize + (New_Y - 1)].GetNumOfElementInGraph(), Nodes[Node_name].GetNumOfElementInGraph()) && (Nodes[((New_X - 1) + 1) * PlotSize + (New_Y - 1)].GetOwner() == Pl_color))			// Case for neighborhood node with coordinates (x+1, y-1) AND neighborhood is exists AND edge doesn't exist yet AND player's color == neighborhood's color
			this->AddEdge(Nodes[((New_X - 1) + 1) * PlotSize + (New_Y - 1)].GetNumOfElementInGraph(), Nodes[Node_name].GetNumOfElementInGraph(), Pl_color);		// Make an edge in a graph between these nodes for cost "Pl_color"
		
	}


	void SwapElOnPlot(const int& name_A, const int& name_B)
	{
		// Part 1 - change nodes on the hexplot
		vector <int> oldXY_A(3), oldXY_B(3);

		oldXY_A[1] = Nodes[name_A].GetX();
		oldXY_A[2] = Nodes[name_A].GetY();

		oldXY_B[1] = Nodes[name_B].GetX();
		oldXY_B[2] = Nodes[name_B].GetY();

		HexVertex temp = Nodes[name_A];				// THIS COULD CAUSE TO INCREASE RUNNING TIME. CHECK!
		Nodes[name_A] = Nodes[name_B];
		Nodes[name_B] = temp;

		Nodes[name_A].SetXYandName(oldXY_A[1], oldXY_A[2], PlotSize);
		Nodes[name_B].SetXYandName(oldXY_B[1], oldXY_B[2], PlotSize);

		// Part 2 - change nodes edges (connections) in graph
		this->DeleteEdges(Nodes[name_A].GetNumOfElementInGraph());
		this->DeleteEdges(Nodes[name_B].GetNumOfElementInGraph());

		this->CreateNewEdges(name_A);
		this->CreateNewEdges(name_B);

	}

	// Returns array with all possible names of moves ( zero-element - number of possible moves)
	vector <int> GetAvailableMoves()
	{
		int k = 1;
		vector <int> AvailableMoves(max_num_of_nodes - active_num_of_nodes + 1);		// Num of elements = mun of empty "nodes" on  hexplot

		for (int i = 1; i <= PlotSize; i++)
			for (int j = 1; j <= PlotSize; j++) {

				if (!Nodes[(i - 1) * PlotSize + j].GetName()) {
					AvailableMoves[k] = (i - 1) * PlotSize + j;
					AvailableMoves[0]++;
					k++;
				}

			}

		return AvailableMoves;
	}

	// Shuffles fully filled HexPlot with 
	void Shuffle(vector <int> Moves, const int& non_changable_el_num, const int& num_of_swaps)
	{
		int a, b;

		for (int i = 1; i <= num_of_swaps; i++)
		{
			if ((a = rand() % Moves[0] + 1) == non_changable_el_num)
				if (a != Moves[0])
					a++;
				else
					a--;

			if ((b = rand() % Moves[0] + 1) == non_changable_el_num)
				if (b != Moves[0])
					b++;
				else
					b--;
						
			this->SwapElOnPlot(Moves[a], Moves[b]);
		}

	}


	void PrintNodeInfo(const int& X, const int& Y)
	{
		int name = (X - 1) * PlotSize + Y;

		cout << "Node name: " << Nodes[name].GetName() << endl;
		cout << "Node X-coord: " << Nodes[name].GetX() << endl;
		cout << "Node Y-coord: " << Nodes[name].GetY() << endl;
		cout << "Node Owner: " << Nodes[name].GetOwner() << endl;
		cout << "Node number in Graph: " << Nodes[name].GetNumOfElementInGraph() << endl;
		cout << "Node connections in Graph: ";
		this->Neighbors(Nodes[name].GetNumOfElementInGraph());
		
	}

	//OVERLOADING OPERATORS FOR CLASS HEXPLOT
	//
	// Overloading << operator - print out a whole Hex playing plot with current positions of the players.
	friend ostream& operator<<(ostream& out, HexPlot& MainPlot)
	{

		int k = 1;			// Total lines in hexplot
		int i = 0;
		string Sp = " ";
		int curr_Y = 1;		// Columns identificator

		for (int i = 0; i <= MainPlot.GetPlotSize(); i++)
			if (i == 0)
				cout << "    ";
			else if (i < 9)
				cout << i << "   ";
			else
				cout << i << "  ";


		while (k != 2 * MainPlot.GetPlotSize() + 1) {

			if (!(k % 2)) {

				for (int j = 0; j <= MainPlot.GetPlotSize(); j++) {

					if (j == 0) {

						cout << curr_Y++;
						if (curr_Y >= 11)
							cout << " ";
						else
							cout << "  ";

					}

					// Prints out node and its owner
					else {

						HexVertex* currNode = MainPlot.GetNodeInfoByXY(i, j);

						if (currNode->GetName() == 0)
							cout << "o";
						else {
							if (currNode->GetOwner() == 1)
								cout << "B";
							else
								cout << "W";
						}

						if (j != MainPlot.GetPlotSize())
							cout << "---";

					}

				}

				k++;
				cout << endl << Sp;
				Sp += " ";

			}
			else {

				for (int j = 0; j <= MainPlot.GetPlotSize(); j++) {

					if ((j != 0) && (i != 0)) {

						if (j != MainPlot.GetPlotSize())
							cout << "\\ / ";
						else
							cout << "\\";

					}
					else
						cout << "   ";

				}

				i++;
				k++;
				cout << endl << Sp;
				Sp += " ";

			}

		}

		return out;

	}


};


class Player
{
private:

	HexPlot PlayerPlot;			// Hexplot of the player with his color nodes only
	int Color;					// What color has a player: "B" - 1 or "W" - 2) 
	string Name;				// Player's name
	int PlotSize;

	bool is_Win;				// Does the player has a winning position

	bool Check_XY(vector <int>& XY_coord, HexPlot& MainPlot)
	{
		if (MainPlot.GetNodeInfoByXY(XY_coord[0], XY_coord[1])->GetName())
			return 0;		// Such node has already exists
		else
			return 1;		// It's a new node
	}


public:

	//CONSTRUCTORS FOR CLASS PLAYER
	//
	// Constructor by default
	Player()
	{
		PlotSize = 2;
		PlayerPlot.SetPlotSize(2);
		Color = 1;
		is_Win = 0;
		Name = "Player";
	}

	// Main constructor - take Player color and hexplot size
	Player(string Pl_name, const int& Pl_Color, const int& PSize)
	{
		PlotSize = PSize;
		PlayerPlot.SetPlotSize(PlotSize);
		Color = Pl_Color;
		is_Win = 0;
		Name = Pl_name;
	}

	int GetPlayerColor() { return Color; }

	string GetPlayerName() { return Name; }

	HexPlot& GetPlayerPlot() { return PlayerPlot; }
	

	// Makes a move for a Player and checks whether the move has been a winnig move (storing this result in "is_Win" var)
	void MakeMove(const int& New_X, const int& New_Y)
	{
		if (PlayerPlot.NewMove(New_X, New_Y, Color) == 1) {

			//Check for winning if there are more than size moves
			is_Win = PlayerPlot.isWin(New_X, New_Y, Color);
		}
	}


	// Main method for AI. It generates optimal X-Y for the next move.
	vector <int> GenerateMove(HexPlot& MainPlot)
	{
		const int NUM_OF_GAMES = 1000;		

		int
			currX,									// X for current analyzing move 
			currY,									// Y for current analyzing move
			num_of_wins,							// How many wins we have for currX-currY move
			shuffle_coefficient,					// How many times elements on theoretical filled hexplot should be swapped. The bigger - the better, but the bigger - more slow move would be generated
			findingX,								// Auxilary X for finding node with Player2 (White)
			findingY,								// Auxilary Y for finding node with Player2 (White)
			j;										// Counter

		vector <int> 
			PossibleMoves = MainPlot.GetAvailableMoves(),			// vector with size in [0] and all names of free dots in its elements			
			BestXY(3);												// Best move among all possible, according to MonteCarlo random generation method. [0] - num of wins according to MC-method, [1] is X for move, [2] is Y for move

		HexPlot tempPlot;
		tempPlot = MainPlot;

		tempPlot.MakeRandomDistr();									// Generates first random filled hexplot according to moves have done earlier

		//Determine optimal shuffling coefficient for current position on hexplot depending on the amount of possible moves on current hexplot
		if (PossibleMoves[0] >= 70)
			shuffle_coefficient = 70;			
		else if ( PossibleMoves[0] < 15 )
			shuffle_coefficient = 15; 			
		else
			shuffle_coefficient = 50; 			
				
		// Start to evaluate all possible moves for Player2
		for (int c = 1; c <= PossibleMoves[0]; c++)					// Cycle from 1 to number of possible moves
		{
			
			currY = (PossibleMoves[c] % PlotSize) ? (PossibleMoves[c] % PlotSize) : PlotSize;							// Convert move node name to its Y
			currX = (currY == PlotSize) ? (PossibleMoves[c] / PlotSize) : (PossibleMoves[c] / PlotSize + 1);			// Convert move node name to its X
			
			// Determine whether current position for theoretical move is owned by 2nd player (White). If no - try to find first White on hexplot and make a swap between this position and current theoretical position
			if (tempPlot.GetNodeInfoByXY(currX, currY)->GetOwner() != 2) {
				
				j = 0;

				do {

					j++;
				
					findingY = (PossibleMoves[j] % PlotSize) ? (PossibleMoves[j] % PlotSize) : PlotSize;							// Convert move node name to its Y
					findingX = (findingY == PlotSize) ? (PossibleMoves[j] / PlotSize) : (PossibleMoves[j] / PlotSize + 1);			// Convert move node name to its X						

				} while (tempPlot.GetNodeInfoByXY(findingX, findingY)->GetOwner() != 2);

				tempPlot.SwapElOnPlot(PossibleMoves[c], PossibleMoves[j]);

			}

			num_of_wins = 0;

			// Start to play "NUM_OF_GAMES" times with currX-currY first move
			for (int i = 1; i <= NUM_OF_GAMES; i++) {
								
				tempPlot.Shuffle(PossibleMoves, c, shuffle_coefficient);
				
				if (tempPlot.WhoWins() == 2)
					num_of_wins++;

			}

			// If current X-Y move leads to better position in game - BestXY = current X-Y
			if (num_of_wins > BestXY[0]) {
				BestXY[0] = num_of_wins;
				BestXY[1] = currX;
				BestXY[2] = currY;
			}
			
		}

		return BestXY;

	}
	

	// Generate just random move. Needs for testing. DON'T USE IN FINAL PROGRAM
	vector <int> GenerateDummyMove(HexPlot& MainPlot)
	{
		srand(time(0));

		vector <int> XY_coord(3);

		while (true) {
			XY_coord[0] = 0;
			XY_coord[1] = rand() % PlotSize + 1;
			XY_coord[2] = rand() % PlotSize + 1;

			if (Check_XY(XY_coord, MainPlot))
				return XY_coord;
		}

		return XY_coord;
	}


	int CheckForWin() { return is_Win; }

	// Method that takes the particular number of move of the player (= num of vertex in graph) and try to find this node in graph for providing X and Y of this move
	vector<int> GetMoveInfo(HexPlot& MainPlot, const int& num_in_graph)
	{
		vector <int> Str(2);
		
		for (int i = 1; i <= MainPlot.GetPlotSize(); i++)
			for (int j = 1; i <= MainPlot.GetPlotSize(); j++)

				if ( MainPlot.GetNodeInfoByXY(i, j)->GetName() && (MainPlot.GetNodeInfoByXY(i, j)->GetOwner() == Color) ) {			// If (i,j) node exists AND its position in Graph corresonds to Player's moves, THAN...

					if ( ( (Color % 2) && (MainPlot.GetNodeInfoByXY(i, j)->GetNumOfElementInGraph() / 2 + 1 == num_in_graph) )		// Condition for Player1
						||
						( !(Color % 2) && (MainPlot.GetNodeInfoByXY(i, j)->GetNumOfElementInGraph() / 2 == num_in_graph) )			// Condition for Player2
						) 
					{
						Str[0] = MainPlot.GetNodeInfoByXY(i, j)->GetX();
						Str[1] = MainPlot.GetNodeInfoByXY(i, j)->GetY();

						return Str;
					}
																
				}						
	}


	// Print-out info about particular node
	void ShowNodes()
	{
		for (int i = 1; i <= PlayerPlot.GetNumOfActiveNodes(); i++)
			//for (int j = 1; j <=PlayerPlot.GetNumOfActiveNodes(); j++)
			PlayerPlot.Neighbors(i);
	}

};



int main(void)
{
	srand(time(0));
		
	int size = 0,						// Hexplot size
		winner = 0,							// Who wins: 0 - nobody, 1 - Player1, 2 - Player2		
		Move_X,							// X-coordinate of Player1 moving
		Move_Y,							// Y-coordinate of Player1 moving
		turn = 1;						// Who should make a move now (1 - Player1, 2 - Player2)

	string PlayerName;
	vector <int> XY_coord(3);			// Need for storing automatically generated X-Y coordinates for computer's move: XY_coord[1] is X, XY_coord[2] is Y
	
	do {
		cout << "Enter the size of hexplot for gaming = ";
		cin >> size;
		if ((size < 2) || (size > 11))
			cout << "Hexplot should be within 2 and 11. Please re-enter the size" << endl;
	} while ((size < 2) || (size > 11));


	cout << "Enter your name (for Black): ";
	cin >> PlayerName;

	// Creation and initializing main hexplot
	HexPlot MainPlot;
	MainPlot.SetPlotSize(size);

	// Create 2 players
	Player Player1(PlayerName, 1, size), Player2("Compuster", 2, size);

	cout << MainPlot << endl << endl;

	cout << "GAME STARTS!" << endl;
	cin.get();

	
	// GENERAL GAME CYCLE
	while (!winner) {			// Until somebody wins - do this cycle

		system("cls");

		cout << endl;
		cout << "Previous moves:" << endl << endl;
		cout << "   " << Player1.GetPlayerName() << " (Black, N-to-S)" << "\t\t|\t" << Player2.GetPlayerName() << " (White, E-to-W)" << endl << endl;
		if (MainPlot.GetNumOfActiveNodes() > 1)
			for (int i = 1; i <= MainPlot.GetNumOfActiveNodes() / 2; i++) {
				XY_coord = Player1.GetMoveInfo(MainPlot, i);
				cout << i << ". (" << XY_coord[0] << " ; " << XY_coord[1] << ")" << "\t\t\t|\t\t";
				XY_coord = Player2.GetMoveInfo(MainPlot, i);
				cout << "  (" << XY_coord[0] << " ; " << XY_coord[1] << ")" << endl;
				cout << "------------------------------------------------------------------------------------" << endl;
			}

		cout << endl << "Main hexplot:" << endl << endl;
		cout << MainPlot << endl << endl;
		
		if ((turn % 2)) {		// Player's 1 turn to make a move

			do {

				cout << "Now player " << Player1.GetPlayerName() << " (Black) make a move. Enter X-Y coordinates of a node you wish to occupy:" << endl;
				cout << "X = ";
				cin >> Move_X;
				cout << "Y = ";
				cin >> Move_Y;
												
			} while (!MainPlot.NewMove(Move_X, Move_Y, Player1.GetPlayerColor()));			// Make a move and check whether it's legal or not. If not - will make a move one more time

			//Check is it a winnig move or not
			winner = MainPlot.WhoWins();
									
			turn = 2;			// Take the turn to another player
			
		}
		else {			// Computer's turn to make a move

			cout << "Now player " << Player2.GetPlayerName() << " make a move." << endl;			
			cout << "Thinking... Please, wait (no more than 2 minutes)" << endl;

			XY_coord = Player2.GenerateMove(MainPlot); 			
			MainPlot.NewMove(XY_coord[1], XY_coord[2], Player2.GetPlayerColor());
			cout << "Ready!" << endl;

			//Check is it a winnig move or not
			winner = MainPlot.WhoWins();

			turn = 1;		// Take the turn to Player1

		}

	}


	// Finishing screen
	system("cls");

	cout << endl;
	cout << "Previous moves:" << endl << endl;
	cout << "   " << Player1.GetPlayerName() << " (Black, N-to-S)" << "\t\t|\t" << Player2.GetPlayerName() << " (White, E-to-W)" << endl << endl;
	
	// Final all moves table
	int r;	// Temp var

	(MainPlot.GetNumOfActiveNodes() % 2) ? r = MainPlot.GetNumOfActiveNodes() + 1 : r = MainPlot.GetNumOfActiveNodes();

	for (int i = 1; i <= r / 2; i++) {
		XY_coord = Player1.GetMoveInfo(MainPlot, i);
		cout << i << ". (" << XY_coord[0] << " ; " << XY_coord[1] << ")" << "\t\t\t|\t\t";
		if ((winner == 1) && (i == r / 2))
			cout << "  ----" << endl;
		else
		{
			XY_coord = Player2.GetMoveInfo(MainPlot, i);
			cout << "  (" << XY_coord[0] << "; " << XY_coord[1] << ")" << endl;
		}

		cout << "------------------------------------------------------------------------------------" << endl;
	}
	
	// Final hexplot configuration
	cout << endl << "Final hexplot configuration:" << endl << endl;
	cout << MainPlot << endl << endl;

	if ( winner == 1 )
		cout << "WINNER IS " << Player1.GetPlayerName() << "! CONGRATULATIONS!!!" << endl;
	else
		cout << "WINNER IS " << Player2.GetPlayerName() << "! CONGRATULATIONS!!!" << endl;

	cout << "Thank you for the game!" << endl << endl;

	//Technical pause
	cin >> size;

	// END BLOCK FOR A REAL GAME
		
	return 0;

}