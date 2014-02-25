//
//  main.cpp
//  clubs
//
//  Oliver Ernst
//  Released under GPL v2
//
//  Clustering Using Binary Splitting (CLUBS) algorithm implemented in C++
//  Credit to the respective authors
//  Written for TopCoder marathon match; TrainingData is attributed to the respective owners
//

#include <iostream>
#include <string>
#include <vector>
#include <sstream> // string manip
#include <math.h> // math
#include <stdlib.h> // for error
#include <fstream> // Reading data

#include <stdio.h>


using namespace std;

/* Function to check if element x is in vector v
    Input : v = vector
            x = element
    Output: 1 if yes, 0 if no
 */
template<typename T>
int in_vec(vector<T> &v, T x)
{
    if(find(v.begin(), v.end(), x) != v.end())
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

/* Function to find the maximum element in a vector, excluding some indices
    Input: v = vector of double to check
            taken_indices = vector of integers of indexes to exclude
    Output: integer of maximum element in v
 */
int max_vec(vector<double> &v, vector<int> &taken_indices)
{
    // Check that there are enough indices available for this to make sense!
    if (v.size() <= taken_indices.size()) {return -1;}
    
    int max=0;
    for (int k=0; k<v.size(); k++)
    {
        if (!in_vec(taken_indices,k) && v[k] > v[max]) {max = k;}
    }
    return max;
}

/* Function to find the minimum element in a vector, excluding some indices
 Input: v = vector of double to check
 taken_indices = vector of integers of indexes to exclude
 Output: integer of minimum element in v
 */
int min_vec(vector<double> &v, vector<int> &taken_indices)
{
    // Check that there are enough indices available for this to make sense!
    if (v.size() <= taken_indices.size()) {return -1;}
    
    int max=0;
    for (int k=0; k<v.size(); k++)
    {
        if (!in_vec(taken_indices,k) && v[k] < v[max]) {max = k;}
    }
    return max;
}

/* Point structure for sorting
 */
struct point {
    int index;
    double coord;
};

struct sort_by_coord {
    bool operator()(point const &a, point const &b) {
        return a.coord < b.coord;
    }
};

/* Sort a vector of cluster indices by their coordinates in S
 Input: C = cluster of indices to original data points
 S = those original, unordered data points (we are looking up S[C[i]]
 */
void sort_vector(vector<int> &C, vector<double> &S) {
    vector<point> Vsort;
    size_t M=C.size();
    
    // Read in as points
    point a;
    for (size_t i=0; i<M; i++)
    {
        a.index = C[i]; a.coord = S[C[i]];
        Vsort.push_back(a);
    }
    
    // Sort
    sort(Vsort.begin(), Vsort.end(), sort_by_coord());
    
    // Read out
    C.clear();
    for (size_t i=0; i<M; i++)
    {
        C.push_back((Vsort[i]).index);
    }
}

/*
 Sanity check for Cp; check that all the elements only appear once...
 */
int sanity_check_Cp(vector<vector<int> > &Cp, size_t size)
{
    // Check size
    size_t ctr=0;
    for (size_t i=0; i<Cp.size(); i++)
    {
        for (size_t j=0; j<Cp[i].size(); j++)
        {
            ctr++;
        }
    }
    if (ctr != size) {return 0;}
    
    // Check elements
    vector<int> illegal_vector;
    
    for (size_t i=0; i<Cp.size(); i++)
    {
        for (size_t j=0; j<Cp[i].size(); j++)
        {
            if (in_vec(illegal_vector,Cp[i][j])) {return 0;}
            else {illegal_vector.push_back(Cp[i][j]);}
        }
    }
    return 1;
}

/* Function to initialize the tree to cluster
    Input: DS_size = data size
            BT = tree structure (empty)
 */
void initializeTree(size_t DS_size, vector<vector<int> >& BT)
{
    vector<int> tmp;
    for (int i=0; i<DS_size; i++)
    {
        tmp.push_back(i);
    }
    BT.push_back(tmp);
}

/*
 Compute SSQ for a cluster C:
 Input: SSQ = output (written over)
 C = cluster
 SxPts = complete set of points (ie original unordered data)
 SyPts = ''
 QxPts = ''
 QyPts = ''
 */
void computeSSQ(vector<double> &SSQ, vector<int> &C, vector<double> &SxPts, vector<double> &SyPts, vector<double> &QxPts, vector<double> &QyPts)
{
    SSQ.clear();
    size_t N = C.size();
    for (int i=0; i<N; i++)
    {
        SSQ.push_back(QxPts[C[i]] - (pow(SxPts[C[i]],2)/N) + QyPts[C[i]] - (pow(SyPts[C[i]],2)/N));
    }
}

/*
 Average an SSQ vector
 Input: SSQ = SSQ vector (actually, any vector)
 */
double computeAverageDeltaSSQ(vector<double> &SSQ)
{
    size_t M = SSQ.size();
    double sum=0.0;
    for (int i=0; i<M; i++)
    {
        sum += SSQ[i];
    }
    sum /= M;
    return sum;
}

/* Compute change in SSQ from splitting
 Input: j = index of splitting
 Cdim_ordered = set of point indices to original data SxPts and SyPts, but ordered along a particular dimension for a particular cluster
 Sdim = complete set of points (ie original unordered data) along i dimension
 SOdim = '' along other dimension ''
 */
double computeDeltaSSQ(vector<int> &Cdim1_ord, vector<int> &Cdim2_ord, vector<double> &Sdim, vector<double> &SOdim)
{
    // Cdim_ordered = ordered indices of points in Sdim
    // Sdim = SPts along dim (not ordered)
    // SOdim = pts along other dim (not ordered)
    
    // j = split pos is included as the element of the first set
    // j=0 => C1 has 0th element, C2 has rest
    // j=1 => C1 has 0,1 elements, C2 has rest
    // ...
    // j=M-2 => C1 has 0,...M-2 elements, C2 has only (M-1)st final
    
    size_t N1 = Cdim1_ord.size(); // j=0 => N1 = 1
    size_t N2 = Cdim2_ord.size(); // j = C.size()-2 => N2 = 1
    
    // Compute the S1 and S2 sums for each cluster
    // Along the ordered dim
    double S1dim=0.0,S2dim=0.0;
    for (size_t i=0; i<N1; i++) {S1dim+=Sdim[Cdim1_ord[i]];}
    for (size_t i=0; i<N2; i++) {S2dim+=Sdim[Cdim2_ord[i]];}
    
    // Along the other dim
    double S1Odim=0.0,S2Odim=0.0;
    for (size_t i=0; i<N1; i++) {S1Odim+=SOdim[Cdim1_ord[i]];}
    for (size_t i=0; i<N2; i++) {S2Odim+=SOdim[Cdim2_ord[i]];}
    
    return (1.0*N1*N2/(1.0*N1+1.0*N2))*(pow((S1dim/N1)-(S2dim/N2),2) + pow((S1Odim/N1)-(S2Odim/N2),2));
}

/* Top down splitting
 */
double topDownSplitting(vector<vector<int> > &BT, vector<vector<int> > &Cp, vector<double> &SxPts, vector<double> &SyPts, vector<double> &QxPts, vector<double> &QyPts)
{
    // Priority queue
    vector<int> PQ;
    // Finished boolean
    bool FINISHED=false;
    // delta SSQ
    double ave_dssq;
    
    // Add the root to the PQ
    PQ.push_back(0);
    
    // initialize clusters
    Cp.push_back(BT[0]);
    
    // Compute the average SSQ
    vector<double> SSQ;
    computeSSQ(SSQ,Cp[0],SxPts,SyPts,QxPts,QyPts);
    ave_dssq = computeAverageDeltaSSQ(SSQ);
    vector<double> ave_dssq_vec;
    ave_dssq_vec.push_back(ave_dssq);
    double ave_dssq_fixed = ave_dssq;
    
    // Sorted cluster indices
    vector<int> Cx_ordered, Cy_ordered;
    
    size_t cluster_index = 0; // index of cluster with highest ssq; avoids sorting
    // For the maximum delta ssq
    double max_dssq=0.0, split1=0.0, split2=0.0;
    size_t j_max, i_max;
    // For the splitting, for checking clusters
    vector<int> Cx1, Cx2, Cy1, Cy2;
    // For the splitting, the new clusters
    vector<int> Cp1,Cp2;
    
    // size_t it_ctr=0;
    while (!FINISHED)
    {
        // cout << "-----------------" << endl;
        // cout << "Iteration: " << it_ctr++ << endl;
        
        // Find the cluster with the highest ssq
        cluster_index=0;
        for (size_t i=1; i<ave_dssq_vec.size(); i++)
        {
            if (ave_dssq_vec[i]>ave_dssq_vec[cluster_index]) {cluster_index = i;}
        }
        
        // Sort this cluster
        Cx_ordered = Cp[cluster_index]; Cy_ordered = Cp[cluster_index];
        sort_vector(Cx_ordered,SxPts);
        sort_vector(Cy_ordered,SyPts);
        
        // Find the maximum delta ssq
        max_dssq=0.0;
        Cx1.clear(); Cx2.clear(); Cy1.clear(); Cy2.clear();
        // Start with all in Cp2, none in Cp1, then transfer over on at a time
        Cx2=Cx_ordered; Cy2=Cy_ordered;
        for (size_t j=0; j<Cp[cluster_index].size()-1; j++)
        {
            // Increase the partition index j
            Cx1.push_back(Cx2[0]); Cy1.push_back(Cy2[0]);
            Cx2.erase(Cx2.begin()); Cy2.erase(Cy2.begin());
            
            split1=computeDeltaSSQ(Cx1,Cx2,SxPts,SyPts);
            split2=computeDeltaSSQ(Cy1,Cy2,SyPts,SxPts);
            
            if (split1 > max_dssq)
            {
                max_dssq = split1;
                i_max = 0; // x
                j_max = j;
            }
            if (split2 > max_dssq)
            {
                max_dssq = split2;
                i_max = 1; // y
                j_max = j;
            }
        }
        
        // Weighting
        max_dssq = pow(max_dssq,0.94);
        
        // Check for splitting condition
        //if (max_dssq > ave_dssq_vec[cluster_index])
        if (max_dssq > ave_dssq_fixed)
        {
            // Do the split
            Cp1.clear();
            Cp2.clear();
            // Find the two splitting clusters
            for (size_t j=0; j<Cx_ordered.size(); j++)
            {
                if (j<=j_max) {Cp1.push_back((i_max<1)*Cx_ordered[j]+(i_max>0)*Cy_ordered[j]);}
                else if (j>j_max) {Cp2.push_back((i_max<1)*Cx_ordered[j]+(i_max>0)*Cy_ordered[j]);}
            }
            
            // Update the Cp
            // Erase the old element
            Cp.erase(Cp.begin()+cluster_index);
            ave_dssq_vec.erase(ave_dssq_vec.begin()+cluster_index);
            // The new stuff
            Cp.push_back(Cp1);
            Cp.push_back(Cp2);
            computeSSQ(SSQ,Cp1,SxPts,SyPts,QxPts,QyPts);
            ave_dssq = computeAverageDeltaSSQ(SSQ);
            ave_dssq_vec.push_back(ave_dssq);
            computeSSQ(SSQ,Cp2,SxPts,SyPts,QxPts,QyPts);
            ave_dssq = computeAverageDeltaSSQ(SSQ);
            ave_dssq_vec.push_back(ave_dssq);
        }
        else
        {
            ave_dssq_vec[cluster_index] = 0;
            // FINISHED = true;
        }
        
        // Check if dssq_vec is all zero
        FINISHED=true;
        for (size_t i=0; i<ave_dssq_vec.size(); i++)
        {
            if (ave_dssq_vec[i] > 1.0e-10) {FINISHED=false;}
        }
    }
    
    return ave_dssq_fixed;
}



void bottomUpMerging(vector<vector<int> > &Cp, double ave_dssq_fixed, vector<double> &SxPts, vector<double> &SyPts)
{
    double minInc=-1;
    size_t merge1, merge2;
    
    double ssq_inc;
    
    // Consider all possible pairs for a merge
    for (size_t i=0; i<Cp.size(); i++)
    {
        for (size_t j=0; j<Cp.size(); j++)
        {
            if (i != j)
            {
                ssq_inc=computeDeltaSSQ(Cp[i],Cp[j],SxPts,SyPts);
                if (ssq_inc < minInc || minInc < 0)
                {
                    minInc = ssq_inc;
                    merge1 = i;
                    merge2 = j;
                }
            }
        }
    }
    
    // for the merger:
    vector<int> tmp;
    size_t swap;
    
    while (minInc < ave_dssq_fixed)
    {
        // Do the merger
        // Swap so that merge1 < merge2 (index size)
        if (merge1 > merge2) {swap = merge1; merge1 = merge2; merge2 = swap;}
        tmp.clear();
        tmp.reserve( Cp[merge1].size() + Cp[merge2].size() ); // preallocate memory
        tmp.insert( tmp.end(), Cp[merge1].begin(), Cp[merge1].end() );
        tmp.insert( tmp.end(), Cp[merge2].begin(), Cp[merge2].end() );
        // Delete the bigger index first (merge2)
        Cp.erase(Cp.begin()+merge2);
        // And the smaller
        Cp.erase(Cp.begin()+merge1);
        // Add tmp
        Cp.push_back(tmp);
        
        // Get the new best pair
        minInc=-1;
        for (size_t i=0; i<Cp.size(); i++)
        {
            for (size_t j=0; j<Cp.size(); j++)
            {
                if (i != j)
                {
                    ssq_inc=computeDeltaSSQ(Cp[i],Cp[j],SxPts,SyPts);
                    if (ssq_inc < minInc || minInc < 0)
                    {
                        minInc = ssq_inc;
                        merge1 = i;
                        merge2 = j;
                    }
                }
            }
        }
    }
}



class OctaveClassifier {
    string *trainingData, *testData;
public:
    int* classify(string *trainingData, string *testData)
    {
        // Timer
        clock_t start;
        double duration;
        start = std::clock();
        
        // Parse the data into a data vector
        vector<double> SxPts,SyPts;
        string id_s, dataid_s, x_s, y_s, c_s, cat_s;
        int id;
        
        size_t DS_size = 0;
        while (!trainingData[DS_size].empty())
        {
            istringstream ss(trainingData[DS_size++]);
            getline(ss,id_s,',');
            getline(ss,dataid_s,',');
            getline(ss,x_s,',');
            getline(ss,y_s,',');
            getline(ss,c_s,',');
            getline(ss,cat_s,',');
            id = atoi(id_s.c_str());
            if (id == 6733)
            {
                SxPts.push_back(atof(x_s.c_str()));
                SyPts.push_back(atof(y_s.c_str()));
            }
        }
        
        // Check that SxPts is not empty!!
        if (SxPts.size() < 1) {cerr << "Empty!" << endl; exit(1);}
        
        // Precompute Qx, Qy for the pts
        vector<double> QxPts,QyPts;
        for (int i=0; i<DS_size; i++)
        {
            QxPts.push_back(SxPts[i]*SxPts[i]);
            QyPts.push_back(SyPts[i]*SyPts[i]);
        }
        
        // Initialize the tree
        vector<vector<int> > BT;
        initializeTree(DS_size,BT);
        
        // Top down splitting
        vector<vector<int> > Cp;
        double ave_dssq_fixed = topDownSplitting(BT,Cp,SxPts,SyPts,QxPts,QyPts);
        
        // Bottom up merging
        if (Cp.size() > 1)
        {
            bottomUpMerging(Cp,ave_dssq_fixed,SxPts,SyPts);
        }
        
        // Identify clusters:
        // Compute mid points for all clusters
        // => max x and max y clusters
        // Compute midpoint of max x and max y cluster
        // Compute distance of other clusters to this midpoint
        // => other two clusters
        // all others = aberrant
        
        int *point_cat = new int[DS_size];
        
        // Compute the midpoints
        vector<double> xmid,ymid;
        double sumx=0.0,sumy=0.0;
        for (size_t i=0; i<Cp.size(); i++)
        {
            sumx=0.0;sumy=0.0;
            for (size_t j=0; j<Cp[i].size(); j++)
            {
                sumx += SxPts[Cp[i][j]];
                sumy += SyPts[Cp[i][j]];
            }
            sumx/=Cp[i].size(); sumy/=Cp[i].size();
            xmid.push_back(sumx);
            ymid.push_back(sumy);
        }
        
        vector<int> taken_indices;
        
        // Find the maximum x and y
        int indy = max_vec(ymid,taken_indices);
        taken_indices.push_back(indy);
        int indx = max_vec(xmid,taken_indices);
        taken_indices.push_back(indx);
        
        // The midpoint of these
        double midx = 0.5*(xmid[indx] + xmid[indy]);
        double midy = 0.5*(ymid[indx] + ymid[indy]);
        
        // The distance of all the vectors to the midpoint
        vector<double> dist;
        for (size_t i=0; i<Cp.size(); i++)
        {
            dist.push_back(pow(xmid[i]-midx,2)+pow(ymid[i]-midy,2));
        }
        
        // The undetermined is worth more, so try that first
        // First create a box with corners at the x and y clusters; this thing must lie outside that box!
        //vector<double> dist_box;
        //for (size_t i=0; i<Cp.size(); i++)
        int indund = max_vec(dist,taken_indices);
        taken_indices.push_back(indund);
        
        // Both
        int indboth = min_vec(dist,taken_indices);
        
        // All the other are abberant

        // Fill out the int vector
        if (indy >= 0)
        {
            for (size_t i=0; i<Cp[indy].size(); i++) {point_cat[Cp[indy][i]] = 1;}
        }
        if (indx >= 0)
        {
            for (size_t i=0; i<Cp[indx].size(); i++) {point_cat[Cp[indx][i]] = 2;}
        }
        if (indund >= 0)
        {
            for (size_t i=0; i<Cp[indund].size(); i++) {point_cat[Cp[indund][i]] = 3;}
        }
        if (indboth >= 0)
        {
            for (size_t i=0; i<Cp[indboth].size(); i++) {point_cat[Cp[indboth][i]] = 4;}
        }
        for (size_t i=0; i<Cp.size(); i++)
        {
            if (i != indy && i != indx && i != indund && i != indboth)
            {
                for (size_t j=0; j<Cp[i].size(); j++)
                {
                    point_cat[Cp[i][j]] = 5;
                }
            }
        }
        
        // Closing
        
        if (!sanity_check_Cp(Cp,BT[0].size())) {cerr << "ERROR! Sanity check failed; exiting..." << endl; exit(1);}
        
        cout << "----------" << endl;
        
        cout << "Total number of clusters: " << Cp.size() << endl;
        for (size_t k=0; k<Cp.size(); k++)
        {
            cout << "---" << endl;
            cout << "Cluster no. :: " << k << endl;
            cout << "Elements :: ";
            for (size_t l=0; l<Cp[k].size(); l++)
            {
                cout << Cp[k][l] << "  ";
            }
            cout << endl;
        }
        
        // If desired, write out
        /*
        cout << "Writing out..." << endl;
        ofstream f("/Users/oernst/TopCoder/clubs/data.txt");
        for (size_t i1 = 0; i1<Cp.size(); i1 ++)
        {
            for (size_t i2 =0; i2<Cp[i1].size(); i2++)
            {
                f << SxPts[Cp[i1][i2]] << ' ';
            }
            f << '\n';
            for (size_t i2 =0; i2<Cp[i1].size(); i2++)
            {
                f << SyPts[Cp[i1][i2]] << ' ';
            }
            f << '\n';
        }
        
        f.close();
         */
        
        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        cout << "Exectution time: " << duration << endl;
        
        return point_cat;
    }
};


//////////////////////////////////////////
// main function
//////////////////////////////////////////

int main(int argc, const char * argv[])
{
    //////////////////////////////////////////
    // Read in some data
    //////////////////////////////////////////

    int MAX = 500;
    
    string *trainingData,*testData;
    trainingData = new string[MAX];
    testData = new string[1];
    
    // variables for reading in lines
    string line,id_s,dataid_s, x_s, y_s, c_s, cat_s;
    int id;
    int line_ctr=0;
    // For storing the category
    int cat_i;
    int *catdat = new int[MAX];
    
    // Read in some data!
    ifstream file1( "TrainingData.csv" );
    while (getline(file1,line))
    {
        // Read in the line
        istringstream ss(line);
        getline(ss,id_s,',');
        getline(ss,dataid_s,',');
        getline(ss,x_s,',');
        getline(ss,y_s,',');
        getline(ss,c_s,',');
        getline(ss,cat_s);
        
        id = atoi(id_s.c_str());
        if (id == 6733) // Read in this id
        {
            // Save the line to pass to the function in proper format
            trainingData[line_ctr] = line;
            
            if (!cat_s.compare(0,1,"F")) {cat_i = 1;}
            else if (!cat_s.compare(0,1,"T")) {cat_i = 2;}
            else if (!cat_s.compare(0,1,"B")) {cat_i = 3;}
            else if (!cat_s.compare(0,1,"U")) {cat_i = 4;}
            else if (!cat_s.compare(0,1,"A")) {cat_i = 5;}
            catdat[line_ctr]=cat_i;
            line_ctr++;
            if (line_ctr > MAX) {break;}
        }
    }
    
    file1.close();
    
    //////////////////////////////////////////
    // Do the clustering
    //////////////////////////////////////////
    
    OctaveClassifier oc;
    int *point_cat = oc.classify(trainingData,testData);
    
    // Score
    double score=0.0;
    for (size_t i=0; i<MAX; i++)
    {
        if (catdat[i]==1)
        {
            if (point_cat[i]==1) {score += 1000;}
            if (point_cat[i]==2) {score += 0;}
            if (point_cat[i]==3) {score += 250;}
            if (point_cat[i]==4) {score += 450;}
            if (point_cat[i]==5) {score += 500;}
        }
        if (catdat[i]==2)
        {
            if (point_cat[i]==1) {score += 0;}
            if (point_cat[i]==2) {score += 1000;}
            if (point_cat[i]==3) {score += 250;}
            if (point_cat[i]==4) {score += 450;}
            if (point_cat[i]==5) {score += 500;}
        }
        if (catdat[i]==3)
        {
            if (point_cat[i]==1) {score += 350;}
            if (point_cat[i]==2) {score += 350;}
            if (point_cat[i]==3) {score += 1000;}
            if (point_cat[i]==4) {score += 500;}
            if (point_cat[i]==5) {score += 750;}
        }
        if (catdat[i]==4)
        {
            if (point_cat[i]==1) {score += 500;}
            if (point_cat[i]==2) {score += 500;}
            if (point_cat[i]==3) {score += 600;}
            if (point_cat[i]==4) {score += 1000;}
            if (point_cat[i]==5) {score += 700;}
        }
        if (catdat[i]==5)
        {
            if (point_cat[i]==1) {score += 500;}
            if (point_cat[i]==2) {score += 500;}
            if (point_cat[i]==3) {score += 750;}
            if (point_cat[i]==4) {score += 700;}
            if (point_cat[i]==5) {score += 1000;}
        }
    }
    
    cout << "Final score: " << score*1000/(line_ctr-1) << endl;
    
    delete [] trainingData; delete [] testData; //freed memory
    trainingData = NULL; testData=NULL; //pointed dangling ptr to NULL
    
    return 0;
}

