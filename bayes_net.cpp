#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <math.h> 
#include <iomanip>

// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

// Our graph consists of a list of nodes where each node is represented as follows:
class Graph_Node{

private:
	string Node_Name;  // Variable name 
	vector<int> Children; // Children of a particular node - these are index of nodes in graph.
	vector<int> Parents; // Parents of a particular node- note these are names of parents
	int nvalues;  // Number of categories a variable represented by this node can take
	vector<string> values; // Categories of possible values
	vector<float> CPT; // conditional probability table as a 1-d array . Look for BIF format to understand its meaning
	pair<int, int> position; //cooardinates of nodes in a 2D - graph
public:
	// Constructor- a node is initialised with its name and its categories
    Graph_Node(string name,int n,vector<string> vals, pair<int, int> pos)
	{
		Node_Name=name;
	
		nvalues=n;
		values=vals;
		position = pos;
	}

	string get_name()
	{
		return Node_Name;
	}

	vector<int> get_children()
	{
		return Children;
	}

	vector<int> get_Parents()
	{
		return Parents;
	}

	vector<float> get_CPT()
	{
		return CPT;
	}

	int get_nvalues()
	{
		return nvalues;
	}

	vector<string> get_values()
	{
		return values;
	}

	pair<int, int> get_position()
	{
		return position;
	}  

	void set_CPT(vector<float> new_CPT)
	{
		CPT.clear();
		CPT=new_CPT;
	}

    void set_Parents(vector<int> Parent_Nodes)
    {
        Parents.clear();
        Parents=Parent_Nodes;
    }

	void set_position(pair<int, int> new_position)
	{
		position = new_position;
	}

	int add_Parent(int p)
	{
        for(int i=0;i<Parents.size();i++)
        {
            if(Parents[i] == p)
                return 0;
        }
        Parents.push_back(p);
        return 1;
	}

    // add another node in a graph as a child of this node
    int add_child(int new_child_index )
    {
        for(int i=0;i<Children.size();i++)
        {
            if(Children[i]==new_child_index)
                return 0;
        }
        Children.push_back(new_child_index);
        return 1;
    }

	int indexOfValue(string v)
	{
		for(int i = 0;i < values.size();i++)
		{
			if(values[i].compare(v) == 0)
				return i;
		}
		return -1;
	}


};


 // The whole network represted as a list of nodes
class network{

	vector<Graph_Node> Pres_Graph;

public:
	int push_back(Graph_Node node)
	{
		Pres_Graph.push_back(node);
		return 0;
	}
    
	int size()
	{
		return Pres_Graph.size();
	}

    // get the index of node with a given name
    int get_index(string val_name)
    {
        for(int i = 0;i < Pres_Graph.size();i++)
        {
            if(Pres_Graph[i].get_name().compare(val_name)==0)
                return i;
        }
        return -1;
    }

    // get the node at nth index
    Graph_Node* get(int n)
    {
        return &Pres_Graph[n];
    }

	void printCPT(int node){
		Graph_Node g = Pres_Graph[node];
		
		vector<float> cpt = g.get_CPT();
		vector<int> parents = g.get_Parents();
		// for(int i = 0;i<parents.size();i++)
		// {
		// 	cout<<"i = "<<i<<" parent = "<<Pres_Graph[parents[i]].get_name()<<endl;
		// 	if(node == parents[i])
		// 		cout <<"Shiteeeeeee!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		// }
		// cout<<"Parents size = "<<parents.size()<<flush<<endl;
		int c = 0;
		int f = 1;
		for(int j = 0;j < parents.size();j++)
		{
			int fac = (Pres_Graph[parents[j]].get_nvalues());
			// cout<<"value of fac = "<<fac<<endl;
			f = f * fac;
		}
		// cout<<"value of f1 = "<<f<<endl;
		int f1 = f;
		// int f = pow(2.0, parents.size());
		int x;

		vector<vector<pair<string, string> > > cptNames(cpt.size());
		// cout<<"cpt size ="<<cpt.size()<<endl;
		// for(int i = 0;i<parents.size();i++)
		// {
		// 	cout<<"i = "<<i<<" parent = "<<Pres_Graph[parents[i]].get_name()<<endl;
		// 	if(node == parents[i])
		// 		cout <<"Shiteeeeeee!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		// }

		for(int i = 0;i < cpt.size();i++)
		{
			f = f1;
			x = c / f;
			// cout<<"i = "<<i<<" value of x = "<<x<<endl;
			cptNames[i].push_back(make_pair(g.get_name() ,g.get_values()[x]));
			// cout<<"Node  = "<<g.get_name()<<endl;

			for(int j = 0;j < parents.size();j++)
			{
				int fac = Pres_Graph[parents[j]].get_nvalues(); 

				f  = f / fac;
				x = (c / f) % (fac);
				Graph_Node p = Pres_Graph[parents[j]];
				// cptNames[i].push_back(p.get_values()[x]);
				cptNames[i].push_back(make_pair(p.get_name() ,p.get_values()[x]));
				// cout<<"Parent  = "<<p.get_name()<<endl;
			}
			// cout<<"i = "<<i<<" Size = "<< cptNames[i].size()<<endl;
			c++;
		}
		for(int i = 0;i < cptNames.size();i++)
		{
			cout<< "P( ";
			int j = 0;
			cout << cptNames[i][j].first << " = "<<cptNames[i][j].second;
			// cout<<cptNames[i].size()<<endl;
			if (cptNames[i].size() > 1)
			{
				cout <<" | ";
				for(j = 1;j < cptNames[i].size() - 1;j++)
					cout << cptNames[i][j].first << " = "<<cptNames[i][j].second <<", ";
				cout << cptNames[i][j].first << " = "<<cptNames[i][j].second;
			}			
			cout<<" ) = "<<cpt[i]<<endl;
			// cout<<cpt[i]<<" ";
		}
		cout<<endl;

	}

	float probab(int node, int value, vector<int> parents)
	{
		Graph_Node g = Pres_Graph[node];
		int c = 0;
		int fac = 1;
		for(int i = parents.size() - 1 ;i >= 0 ;i--)
		{
			c += parents[i] * fac;
			fac *= Pres_Graph[g.get_Parents()[i]].get_nvalues();
		}
		c += value * fac;
		return g.get_CPT()[c];
	}

	int CPTindex(int node, vector<int> data)
	{
		int value = data[node];
		Graph_Node g = Pres_Graph[node];
		vector<int> parents = g.get_Parents();
		vector<int> parentValues;
		for(int i = 0;i < parents.size();i++)
		{
			parentValues.push_back(data[parents[i]]);
		}
		int c = 0;
		int fac = 1;
		for(int i = parentValues.size() - 1 ;i >= 0 ;i--)
		{
			c += parentValues[i] * fac;
			fac *= Pres_Graph[g.get_Parents()[i]].get_nvalues();
		}
		c += value * fac;
		return c;
	}


};

network read_network(const char* networkPath)
{
	network Alarm;
	// vector<Graph_Node> nodeNet;
	string line;
	int find=0;
	ifstream myfile(networkPath);
  	// ifstream myfile("alarm.bif"); 
  	string temp;
  	string name;
  	vector<string> values;
  	
    if (myfile.is_open())
    {
    	while (! myfile.eof() )
    	{
    		stringstream ss;
      		getline (myfile,line);
      		
      		ss.str(line);
     		ss>>temp;
     		// cout<<temp<<endl;
     		
     		if(temp.compare("variable")==0)
     		{
                    
     				ss>>name;
     				getline (myfile,line);
                   
     				stringstream ss2;
     				ss2.str(line);
     				for(int i=0;i<4;i++)
     				{
     					ss2>>temp;	
     				}

     				values.clear();
     				while(temp.compare("};")!=0)
     				{
     					values.push_back(temp);
     					// cout<<temp<<endl;
     					ss2>>temp;

    				}
					// cout<<"First temp = "<<temp<<endl;
					getline (myfile,line);
     				stringstream ss3;
     				ss3.str(line);

					for(int i = 0;i < 4;i++)
					{
						ss3>>temp;
					}
					temp = temp.substr(1, temp.size() - 2);
					int x = atoi(temp.c_str());
					ss3>>temp;
					temp = temp.substr(0, temp.size() - 2);
					int y = atoi(temp.c_str());

     				Graph_Node new_node(name,values.size(),values,make_pair(x,y));
     				int pos = Alarm.push_back(new_node);
					// nodeNet.push_back(new_node);
     				
     		}

     		else if(temp.compare("probability")==0)
     		{
                    
     				ss>>temp;
     				ss>>temp;
     				
                    int child = Alarm.get_index(temp);
                    ss>>temp;
                    // vector<int> parents;
     				while(temp.compare(")")!=0)
     				{
                        int parent = Alarm.get_index(temp);
						
                        Alarm.get(parent)->add_child(child);
     					Alarm.get(child)->add_Parent(parent);
     					ss>>temp;

    				}
    				getline (myfile,line);
     				stringstream ss2;
                    
     				ss2.str(line);
     				ss2>> temp;
                    
     				ss2>> temp;
                    
     				vector<float> curr_CPT;
                    string::size_type sz;
     				while(temp.compare(";")!=0)
     				{
     					curr_CPT.push_back(atof(temp.c_str()));
     					ss2>>temp;
					}
                    Alarm.get(child)->set_CPT(curr_CPT);

     		}
            else
            {
                
            }    		
    	}

		if(find == 1)
    		myfile.close();
  	}
  	
  	return Alarm;
}

vector<vector<int> > readData(network net, const char* dataPath){
	ifstream myfile(dataPath); 
	string line;
  	string temp;
	vector<vector<int> > ans;
  	// string name;
  	// vector<string> values;
  	int count =  0;
    if (myfile.is_open())
    {
    	while (! myfile.eof() )
    	{
    		stringstream ss;
      		getline (myfile,line);
			ans.push_back(vector<int>());
      		ss.str(line);
			int qIndex  = -1;
     		for(int i = 0;i < net.size();i++)
			{	
				ss>>temp;
				// if(count == 1 && i < 5)
				// 	cout<<"i = "<<i<<" "<<temp<<endl;
				int ind = net.get(i)->indexOfValue(temp); 
				if(ind == -1)
					qIndex = i;
				ans[count].push_back(ind);
			}
			ans[count].push_back(qIndex);
			count++;
		}
	}
	return ans;
}

void intitialize_CPT(network &net)
{
	for(int i = 0;i < net.size();i++)
	{
		int n = net.get(i)->get_CPT().size();
		vector<float> cpt(n);
		int nvalues = net.get(i)->get_nvalues();
		int m = n / nvalues;
		// int m = n / 2;
		for(int j = 0;j < m;j++)
		// for(int j = 0;j < n;j++)
		{
			float sum = 0;
			int ind,k;
			for(k = 0;k < nvalues - 1;k++)
			{
				ind = j + m * k;
				
				cpt[ind] = ((float) rand() / (float) RAND_MAX ) * (float)(1 - sum);
				sum += cpt[ind];

				// if(j < m)
				// 	cpt[j] = (1.0);
				
				// else
				// 	cpt[j] = (0.0);
			}
			ind = j + m * k;
			cpt[ind] = 1 - sum;

		}
		net.get(i)->set_CPT(cpt);
	}
}

void initialize_data(network &net, vector<vector<int> > &data)
{
	for(int i = 0;i < data.size();i++)
	{
		for(int j = 0;j < data[i].size();j++)
		{
			int ind = data[i][data[i].size() - 1];
			int nvalues = net.get(ind)->get_nvalues();
			data[i][ind] = rand() % nvalues;
		}
	}
}


int predictData(network &net, vector<vector<int> > &data)
{
	int count = 0;
	for(int i = 0;i < data.size();i++)
	{
		int node = data[i][data[i].size() - 1];
		vector<int> parents = net.get(node)->get_Parents();
		vector<int> parentValues;
		for(int j = 0;j < parents.size();j++)
		{
			parentValues.push_back(data[i][parents[j]]);
		}
		
		int nvalues = net.get(node)->get_nvalues();
		float max = 0;
		int maxInd = 0;
		for(int j = 0;j < nvalues;j++)
		{
			float p = net.probab(node, j, parentValues);
			if( p > max)
			{
				max = p;
				maxInd = j;
			}
		}
		// cout<<"Prediction is "<<maxInd<<endl;
		if(data[i][node] != maxInd)
		{
			data[i][node] = maxInd;
			count++;
		}

	}	
	return count;
	// cout<<"Changes made = "<<count<<endl;
}

void calculateCPT(network &net, vector<vector<int> > &data, float laplace)
{
	vector<vector<float> > counter;
	for(int i = 0 ;i < net.size();i++)
	{
		int s = net.get(i)->get_CPT().size();
		counter.push_back(vector<float>(s));
		for(int j = 0;j < s;j++)
			counter[i][j] = laplace;
	}

	for(int i =  0;i < data.size();i++)
	{
		for(int j = 0 ;j < net.size();j++)
		{
			int ind = net.CPTindex(j, data[i]);
			counter[j][ind] += 1;
		}
	}
	
	for(int i = 0;i < counter.size();i++)
	{
		int nvalue = net.get(i)->get_nvalues();
		int diff = counter[i].size() / nvalue;
		// diff  = 4 nvalue = 3
		for(int j = 0;j < diff;j++)
		{
			float sum = 0;
			for(int k = 0;k < nvalue;k++)
			{
				sum += counter[i][j + k * diff];
			}
			for(int k = 0;k < nvalue;k++)
			{
				counter[i][j + k * diff] = counter[i][j + k * diff] / sum;
			}
		}
		net.get(i)->set_CPT(counter[i]);		
	}


}


void write_network(network &net)
{
	ofstream outfile;
	outfile.open("solved_alarm.bif");
	outfile<<"// Bayesian Network in the Interchange Format"<<endl;
	outfile<<"// Produced by BayesianNetworks package in JavaBayes"<<endl;
	outfile<<"// Output created Sun Nov 02 17:58:15 GMT+00:00 1997"<<endl;
	outfile<<"// Bayesian network "<<endl;

	outfile << "network \"Alarm\" { //";
	outfile << net.size()<<" variables and "<<net.size()<<" probability distributions"<<endl;
	outfile<<"}"<<endl;
	for(int i = 0;i < net.size();i++)
	{
		Graph_Node g = *net.get(i);
		outfile << "variable  "<< g.get_name() <<" {"<<" //"<<g.get_nvalues()<<" values"<<endl;
		outfile << "\ttype discrete["<< g.get_nvalues()<<"] {" ;
		vector<string> values = g.get_values();
		for(int j = 0;j < g.get_nvalues();j++)
		{
			outfile<<"  "<<values[j];
		}
		outfile<<" };"<<endl;

		outfile << "\tproperty \"position = ("<< g.get_position().first<<", "<<g.get_position().second<<")\" ;"<<endl;
		outfile<<"}"<<endl;	
	}
	for(int i = 0;i < net.size();i++)
	{
		Graph_Node g = *net.get(i);
		outfile << "probability (  "<< g.get_name() ;
		vector<int> parents = g.get_Parents();
		for(int j = 0;j < parents.size(); j++)
		{
			Graph_Node p = *net.get(parents[j]);
			outfile<<"  "<<p.get_name();
		}
		outfile<<" ) { //"<<1+ parents.size()<<" variable(s) and ";
		outfile<<g.get_CPT().size()<<" values"<<endl;
		outfile << "\ttable ";
		vector<float> cpt = g.get_CPT();
		for(int j = 0;j< cpt.size(); j++)
		{
			outfile<<fixed<< setprecision(4) <<cpt[j]<<" ";
		}
		outfile<<";"<<endl<<"}"<<endl;
	}
	outfile.close();


}

int main(int charc, char* argv[])
{
	network Alarm;
	srand(time(0));
	Alarm = read_network(argv[1]);
	vector<vector<int> >  data = readData(Alarm, argv[2]);

    // cout << Alarm.size() <<endl;
	// Graph_Node g1 = *Alarm.get(1);
	// cout<<g1.get_name()<<endl;
	// cout<<g1.get_nvalues()<<endl;
	// cout<<g1.get_position().first<<" "<<g1.get_position().second<<endl;
	// Alarm.printCPT(1);

	float laplace = 0.025;
	// float laplace = atof(argv[3]);
	
	intitialize_CPT(Alarm);
	// initialize_data(Alarm, data);
	// Alarm.printCPT(1);
	// calculateCPT(Alarm, data, laplace);
	int changes = predictData(Alarm, data);
	// cout<<"Changes are "<<changes<<endl;
	int count =  0;
	int epsilon = data.size() / 1000;
	while(changes >= epsilon || count < 5)
	{
		calculateCPT(Alarm, data, laplace);
		changes = predictData(Alarm, data);
		// cout<<"Changes are "<<changes<<endl;
		// Alarm.printCPT(1);
		// cout<<endl;
		if(changes <= epsilon)
			count++;
		else
			count = 0;
	}
	
	write_network(Alarm);
	
// Example: to do something
	// cout<<"Perfect! Hurrah! \n";
	
}




