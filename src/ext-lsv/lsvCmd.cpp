#include "base/abc/abc.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"
#include "sat/cnf/cnf.h"
#include <set>
#include <vector>
#include <list>
#include <queue>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <map>
#include <typeinfo>

using namespace std;
extern "C" {
  extern Aig_Man_t * Abc_NtkToDar(Abc_Ntk_t * pNtk,int fExors,int fRegisters);
  extern Cnf_Dat_t *     Cnf_Derive( Aig_Man_t * pAig, int nOutputs );
}
static int Lsv_CommandPrintNodes(Abc_Frame_t* pAbc, int argc, char** argv);
static int Lsv_CommandMFFC(Abc_Frame_t* pAbc, int argc, char** argv);
static int lsv_or_bidec(Abc_Frame_t* pAbc, int argc, char** argv);

void init(Abc_Frame_t* pAbc) {
  Cmd_CommandAdd(pAbc, "LSV", "lsv_print_nodes", Lsv_CommandPrintNodes, 0);
  Cmd_CommandAdd(pAbc, "LSV", "lsv_print_msfc", Lsv_CommandMFFC, 0);
  Cmd_CommandAdd(pAbc, "LSV", "lsv_or_bidec", lsv_or_bidec, 0);

}

void destroy(Abc_Frame_t* pAbc) {}

Abc_FrameInitializer_t frame_initializer = {init, destroy};

struct PackageRegistrationManager {
  PackageRegistrationManager() { Abc_FrameAddInitializer(&frame_initializer); }
} lsvPackageRegistrationManager;


class Graph{
private:
    int num_vertex;
    std::vector< std::list<int> > AdjList;
    int *color,             // 0:white, 1:gray, 2:black
        *predecessor,
        *distance,          // for BFS()
        *discover,          // for DFS()
        *finish;
public:
    Graph():num_vertex(0){};
    Graph(int N):num_vertex(N){
        // initialize Adj List
        AdjList.resize(num_vertex);
    };
    void AddEdgeList(int from, int to);
        
    void DFS(int Start);
    void DFSVisit(int vertex, int &time);

    void CCDFS(int vertex,int num11,std::map < int,string > nodemap);                // ?拍DFS 
                // ?拍BFS, ?抵?頛臬??函??    void SetCollapsing(int vertex);
                   // ?啣predecessor, 靘?閮潛, ??閬?};
void Graph::AddEdgeList(int from, int to){
    AdjList[from].push_back(to);
}
void Graph::DFS(int Start){
    color = new int[num_vertex];           // ?蔭閮擃?蝵?    discover = new int[num_vertex];
    finish = new int[num_vertex];
    predecessor = new int[num_vertex];

    int time = 0;                          // ???? 憒?銝?b)
    for (int i = 0; i < num_vertex; i++) { 
        color[i] = 0;
        discover[i] = 0;
        finish[i] = 0;
        predecessor[i] = -1;
    }

    int i = Start;
    for (int j = 0; j < num_vertex; j++) { // 瑼Ｘ??raph銝剔?vertex?質?鋡急?撠
        if (color[i] == 0) {               // ?史ertex銝?質, ?脰?隞亥府vertex雿韏琿?銋?撠?            DFSVisit(i, time);
        }
        i = j;                             // j??AdjList摰韏圈?銝?? 蝣箔???ertex?質◤????    }
}

void Graph::DFSVisit(int vertex, int &time){   // 銝?行?vertex鋡怎?曇??舐?? 靘輸脣DFSVisit()
    color[vertex] = 1;                         // ?ertex憛??啗
    discover[vertex] = ++time;                 // ?湔vertex?iscover??
    for (std::list<int>::iterator itr = AdjList[vertex].begin();   // for loop?憭芷
         itr != AdjList[vertex].end(); itr++) {                    // ???拇挾
        if (color[*itr] == 0) {                // ?交?撠?質?ertex
            predecessor[*itr] = vertex;        // ?湔?酥redecessor
            DFSVisit(*itr, time);              // 蝡隞亙雿?啁???韏琿?, ?脣?啁?DFSVisit()
        }
    }
    color[vertex] = 2;                         // ?鈞ertex撌脩???????銋???vertex敺? 撠憛?
    finish[vertex] = ++time;                   // 銝行?逆inish??
}
void Graph::SetCollapsing(int current){
    int root;  // root
    for (root = current; predecessor[root] >= 0; root = predecessor[root]);

    while (current != root) {
        int parent = predecessor[current];
        predecessor[current] = root;
        current = parent;
    }
}

void Graph::CCDFS(int vertex = 0,int num11 = 0, std::map <int,string> nodemap = {}){

    DFS(vertex);      // 
    
    for (int i = 0; i< num_vertex; i++){
        SetCollapsing(i);
    }
    

    int num_cc = 0;
    for (int i = 0; i < num_vertex; i++) {
        if (predecessor[i] < 0 && i >= num11) {
            std::cout << "MSFC " << num_cc++ << ": " << nodemap[i] ;
            for (int j = 0; j < num_vertex; j++) {
                if (predecessor[j] == i) {
                    std::cout << ","<<        nodemap[j] ;
                }
            }
            std::cout << std::endl;
        }
    }
}


void Lsv_NtkPrintNodes(Abc_Ntk_t* pNtk) {
  Abc_Obj_t* pObj;
  int i;
  Abc_NtkForEachNode(pNtk, pObj, i) {
    printf("Object Id = %d, name = %s\n", Abc_ObjId(pObj), Abc_ObjName(pObj));
    Abc_Obj_t* pFanin;
    int j;
    Abc_ObjForEachFanin(pObj, pFanin, j) {
      printf("  Fanin-%d: Id = %d, name = %s\n", j, Abc_ObjId(pFanin),
             Abc_ObjName(pFanin));
    }
    if (Abc_NtkHasSop(pNtk)) {
      printf("The SOP of this node:\n%s", (char*)pObj->pData);
    }
  }
}

int Lsv_CommandPrintNodes(Abc_Frame_t* pAbc, int argc, char** argv) {
  Abc_Ntk_t* pNtk = Abc_FrameReadNtk(pAbc);
  int c;
  Extra_UtilGetoptReset();
  while ((c = Extra_UtilGetopt(argc, argv, "h")) != EOF) {
    switch (c) {
      case 'h':
        goto usage;
      default:
        goto usage;
    }
  }
  if (!pNtk) {
    Abc_Print(-1, "Empty network.\n");
    return 1;
  }
  Lsv_NtkPrintNodes(pNtk);
  return 0;

usage:
  Abc_Print(-2, "usage: lsv_print_nodes [-h]\n");
  Abc_Print(-2, "\t        prints the nodes in the network\n");
  Abc_Print(-2, "\t-h    : print the command usage\n");
  return 1;
}
bool cmp1(std::pair<int,int>a,std::pair<int,int>b)
{
    return a.first < b.first;
}
int Lsv_CommandMFFC(Abc_Frame_t* pAbc, int argc, char** argv) {
  Abc_Ntk_t* pNtk = Abc_FrameReadNtk(pAbc);
  int c;
  Extra_UtilGetoptReset();
  while ((c = Extra_UtilGetopt(argc, argv, "h")) != EOF) {
    
  }
  if (!pNtk) {
    Abc_Print(-1, "Empty network.\n");
    return 1;
  }
  //main
  Abc_Obj_t* pObj;
  int i;
  //static int 14 = Abc_NtkNodeNum(pNtk);
  std::vector <int> if_two_fanout;
  std::vector <std::pair<int,int> > big_set;
  std::vector <int> boom;
  std::map<int, std::string> nodemap;
  int num11 = Abc_NtkObjNum(pNtk) - Abc_NtkNodeNum(pNtk) ;
  //std::cout<<num11;
  //Abc_NtkIncrementTravId(pNtk);
  int node_num = Abc_NtkObjNum(pNtk);
  Abc_NtkForEachNode(pNtk, pObj, i) {
    //printf("Object Id = %d, name = %s\n", Abc_ObjId(pObj), Abc_ObjName(pObj));
    nodemap.insert(std::pair<int, std::string>(Abc_ObjId(pObj), Abc_ObjName(pObj)));
    std::pair<int,int> aig_set;
    aig_set.first=(Abc_ObjId(pObj));
    Abc_Obj_t* pFanin;
    int j;
    Abc_ObjForEachFanin(pObj, pFanin, j) {
      //printf("  Fanin-%d: Id = %d, name = %s\n", j, Abc_ObjId(pFanin),
      //       Abc_ObjName(pFanin));
      if (!Abc_ObjIsPi(pFanin)){

        {
          
          if (std::find(if_two_fanout.begin(), if_two_fanout.end(),Abc_ObjId(pFanin))!=if_two_fanout.end())
          {
            //std::cout<<Abc_ObjId(pFanin)<<std::endl;
            //boom.push_back(Abc_ObjId(pFanin));
            for (int i  = 0; i< big_set.size();i++){
              if(big_set[i].second == Abc_ObjId(pFanin)) big_set.erase(big_set.begin()+i);
            }
          }
          
          else
          {
            aig_set.second = (Abc_ObjId(pFanin));
            big_set.push_back(aig_set);
            if_two_fanout.push_back(Abc_ObjId(pFanin));
          }
          
          
          
        }
      }
    }
    //std::vector<std::set<int> >::iterator it_i_2;
    //std::set <int>::iterator it;
    //for (it = aig_set.begin(); it != aig_set.end(); ++it) {
    //  for(it_i_2=big_set.begin(); it_i_2!=big_set.end(); ++it_i_2)
    //  {
    //    if(*it_i.count(*it) ){
    //      std::set <int> concat_set;
    //      std::set_union(*it_i.begin(), *it_i.end(),
    //            *it.begin(), *it.end(),
    //            std::inserter(concat_set, concat_set.begin()));
    //      
    //    }
    //  }
    //  
    //}
    
    
  }
  //graph start
  
  Graph g3(node_num);


  
//   for(int lp = 0; lp<big_set.size(); lp++){
//     set<int>::iterator it1;
//     for (it1 = big_set[lp].begin(); it1 != big_set[lp].end(); it1 ++)
// {
//     cout << *it1 <<" ";
// }
// std::cout<<std::endl;
//   }
    
  
  //boom.clear();
  if_two_fanout.clear();
  std::vector<std::pair<int, int> > vec;
  for( int op = 0; op<big_set.size(); op++){
    //if (big_set[op].size() >1) {
      //std::cout<<*(big_set[op].begin())<<" "<<*(++(big_set[op].begin()))<<std::endl;
      vec.push_back(std::make_pair(big_set[op].first,big_set[op].second));
      vec.push_back(std::make_pair(big_set[op].second,big_set[op].first));
      
    //}
    //if (big_set[op].size() == 3){
    //  g3.AddEdgeList(*(big_set[op].begin())-num11, *(++(big_set[op].begin()))-num11);
    //  g3.AddEdgeList(*(++(big_set[op].begin()))-num11, *(--(big_set[op].end()))-num11);
    //  //std::cout<<*(big_set[op].begin())<<" "<<*(++(big_set[op].begin()))<<std::endl;
    //  //std::cout<<*(++(big_set[op].begin()))<<" "<<*(--(big_set[op].end()))<<std::endl;
    //} 
    //std::cout<<"Next set"<<std::endl;
        
       // if(itr!=(*it_i).end() && p < (*it_i).size()){    
       //   //g3.AddEdgeList(*itr++, *itr); 
       //   itr--;
       // }
  }
  
  sort(vec.begin(), vec.end(), cmp1);
  for (int re = 0; re<vec.size();re++){
    g3.AddEdgeList(vec[re].first,vec[re].second);
    //std::cout<<vec[re].first<< " "<<vec[re].second<<std::endl;
  }
  big_set.clear();
  //big_set.shrink_to_fit();   
  
  //std::cout<<node_num;
  g3.CCDFS(0,num11,nodemap);
  //std::cout<<node_num<<"  "<<num11<<" ";
  return 0;


}

int lsv_or_bidec(Abc_Frame_t* pAbc, int argc, char** argv) {
    Abc_Ntk_t* pNtk = Abc_FrameReadNtk(pAbc);
    
    Abc_Obj_t* pObj;
    Abc_Ntk_t * pNtkOn1;
    sat_solver * pSat;
    Cnf_Dat_t * pCnf;
    int i;
    Abc_NtkForEachPo(pNtk,pObj,i){
      
      pNtkOn1 = Abc_NtkCreateCone( pNtk, Abc_ObjFanin0(pObj), Abc_ObjName(pObj), 0 );
      if ( Abc_ObjFaninC0(pObj) )
          Abc_ObjXorFaninC( Abc_NtkPo(pNtkOn1, 0), 0 );
      Aig_Man_t* pMan = Abc_NtkToDar(pNtkOn1,0,0);
      Aig_Obj_t * pCo = Aig_ManCo(  pMan,0 ) ;//PO 
      Aig_Obj_t * pObj_Ci;
      pCnf = Cnf_Derive( pMan, 1 );
      int PI_var_num [Aig_ManCiNum(pMan)];
      int PO_var_num = pCnf->pVarNums[pCo->Id];
      int z = 0;
      int j;
      int nInitVars = pCnf->nVars;//varshift
      int nCi = Aig_ManCiNum(pMan);
      
      
      pSat = (sat_solver *)Cnf_DataWriteIntoSolver(pCnf,1,0);
      sat_solver_setnvars(pSat, 3*nInitVars+2*nCi);
      lit Lits[3];
      Lits[0] = toLitCond(PO_var_num,0);
      sat_solver_addclause(pSat,Lits,Lits+1);

      Cnf_Dat_t * pCnf1 = Cnf_DataDup(pCnf);
      Cnf_DataLift(pCnf1,nInitVars);
      for (int i = 0; i < pCnf1->nClauses; i++ ){//changable
        int temp = sat_solver_addclause( pSat, pCnf1->pClauses[i], pCnf1->pClauses[i+1] );
        if ( !temp ){
          sat_solver_delete( pSat );
          std::cout<<"addclause fail!"<<std::endl;
          
        }
      }

      int PO_var_num1 = pCnf1->pVarNums[pCo->Id];
      Lits[0] = toLitCond(PO_var_num1,1);
      sat_solver_addclause(pSat,Lits,Lits+1);
  
      Cnf_Dat_t * pCnf2 = Cnf_DataDup(pCnf1);
      Cnf_DataLift(pCnf2,nInitVars);
      for (int i = 0; i < pCnf2->nClauses; i++ ){//changable
        int temp = sat_solver_addclause( pSat, pCnf2->pClauses[i], pCnf2->pClauses[i+1] );
        if ( !temp ){
          sat_solver_delete( pSat );
          std::cout<<"addclause fail!"<<std::endl;
          
        }
      }

      int PO_var_num2 = pCnf2->pVarNums[pCo->Id];
      Lits[0] = toLitCond(PO_var_num2,1);
      sat_solver_addclause(pSat,Lits,Lits+1);


      int alpha  = 3*(nInitVars)+1,
          //let beta = alpha+1 
          id1st, id2nd, id3rd;

      //Aig_Obj_t* pAigCi;
      Aig_ManForEachCi( pMan, pObj_Ci, j ){
        id3rd = pCnf->pVarNums[pObj_Ci->Id];
        id2nd = id3rd - nInitVars;
        id1st = id2nd - nInitVars;

        //X X' alpha
        //Lsv_AddClauseForEq( pSat, id1st, id2nd, alpha++);
        int Cid;
        //assert( iVarA >= 0 && iVarB >= 0 && iVarEn >= 0 );
        lit Lits2 [3];
        Lits2[0] = toLitCond( id1st, 0 );
        Lits2[1] = toLitCond( id2nd, 1 );
        Lits2[2] = toLitCond( alpha, 0 );
        Cid = sat_solver_addclause( pSat, Lits2, Lits2 + 3 );
        assert( Cid );

        Lits2[0] = toLitCond( id1st, 1 );
        Lits2[1] = toLitCond( id2nd, 0 );
        Lits2[2] = toLitCond( alpha, 0 );
        Cid = sat_solver_addclause( pSat, Lits2, Lits2 + 3 );
        assert( Cid );
        alpha++;
        //X X" beta
        //Lsv_AddClauseForEq( pSat, id1st, id3rd, alpha++);  
        Lits2[0] = toLitCond( id1st, 0 );
        Lits2[1] = toLitCond( id3rd, 1 );
        Lits2[2] = toLitCond( alpha, 0 );
        Cid = sat_solver_addclause( pSat, Lits2, Lits2 + 3 );
        assert( Cid );

        Lits2[0] = toLitCond( id1st, 1 );
        Lits2[1] = toLitCond( id3rd, 0 );
        Lits2[2] = toLitCond( alpha, 0 );
        Cid = sat_solver_addclause( pSat, Lits2, Lits2 + 3 );
        assert( Cid );
        alpha++;
      }
      // for ( int m = 0;m< varshift;m++){
      //       Lits[0] = unit_assumption[m].first;
      //       Lits[1] = toLitCond(PI_var_num[m],1);
      //       Lits[2] = toLitCond(PI_var_num[m]+varshift,0);
      //       sat_solver_addclause(pSat,Lits,Lits+3);
      //       Lits[0] = unit_assumption[m].first;
      //       Lits[1] = toLitCond(PI_var_num[m],0);
      //       Lits[2] = toLitCond(PI_var_num[m]+varshift,1);
      //       sat_solver_addclause(pSat,Lits,Lits+3);
      //     }
      //     for (int m = 0;m< varshift;m++){
      //       Lits[0] = unit_assumption[m].second;
      //       Lits[1] = toLitCond(PI_var_num[m],1);
      //       Lits[2] = toLitCond(PI_var_num[m]+2*varshift,0);
      //       sat_solver_addclause(pSat,Lits,Lits+3);
      //       Lits[0] = unit_assumption[m].second;
      //       Lits[1] = toLitCond(PI_var_num[m],0);
      //       Lits[2] = toLitCond(PI_var_num[m]+2*varshift,1);
      //       sat_solver_addclause(pSat,Lits,Lits+3);
      //     }
      //Aig_ManForEachCi( pMan, pObj_Ci,  j ){
      //  PI_var_num[z] = pCnf->pVarNums[pObj_Ci->Id];
      //  z++;
      //}

      lit Lits1 [2*nCi];
      for (int m = 0;m<2*nCi;m++){
        Lits1[j] = toLitCond( (3*nInitVars + m + 1), 1 );
      }





      //int varshift = sat_solver_nvars(pSat);
      //
      //
      //
      //int PO_var_num1 = PO_var_num + varshift;
      //int PO_var_num2 = PO_var_num1 + varshift;
      
      
      
      
      //Cnf_Dat_t * pCnf2 = Cnf_DataDup(pCnf);
      //Cnf_DataLift(pCnf2,varshift);
      //Cnf_DataLift(pCnf2,varshift);
      //Lits[0] = toLitCond(PO_var_num2,1);
      //sat_solver_addclause(pSat,Lits,Lits+1);
      lbool status = l_True;
      //std::vector<std::pair<int,int>> unit_assumption(varshift, std::make_pair(0, 0));//first = alpha ,second = beta;
      for (int k = 1;k < nCi;++k){
        for (int l = 0;l < k; ++l){
          Lits1[ 2*k   ]--; //alphaj=1
          Lits1[(2*l)+1]--; //betak=1
          
          status = sat_solver_solve(pSat,Lits1,Lits1+2*nCi,0,0,0,0);
          Lits1[ 2*k   ]++; //alphaj=1
          Lits1[(2*l)+1]++; //betak=1
          //final_size = sat_solver_final(pSat,)
          if (status == l_False){
            break;
          }
          
          
        }
        if (status == l_False){
            break;
        }
      }

      if(status==l_False){
        //unsat -> print
        std::cout<<"PO "<<Abc_ObjName(pObj)<<" support partition: 1"<<std::endl;

        //now Lits is alpha,beta
        //initialize to all 1
        for(j=0; j<2*nCi; ++j){ Lits1[j] = 1; }

        int *pClauses, *pLit, *pStop;
        int nFinalClause = sat_solver_final(pSat, &pClauses);
        //std::cerr<<"nFinalClauses: "<<nFinalClause<<std::endl;
        //assert at least two 1 in Lits
        

        for(j=0; j<nFinalClause; ++j){
          for ( pLit = &pClauses[j], pStop = &pClauses[j+1]; pLit < pStop; pLit++ ){
            Lits1[Abc_Lit2Var(*pLit)-(nInitVars*3)-1]=0;
            //std::cout<<Abc_Lit2Var(*pLit)<<" ";
          } //std::cout<<std::endl;
        }

        //just try to make XA and XB balance
        int nXA=0, nXB=0;
        for(j=0; j<nCi; ++j){
          if     (Lits1[(2*j)]==0 && Lits1[(2*j)+1]==1) std::cout<<"1"; //alpha=0, beta=1 => XB
          else if(Lits1[(2*j)]==1 && Lits1[(2*j)+1]==0) std::cout<<"2"; //alpha=1, beta=0 => XA
          else if(Lits1[(2*j)]==0 && Lits1[(2*j)+1]==0) std::cout<<"0"; //alpha=0, beta=0 => XC
          else if(nXA>nXB){ std::cout<<"1"; ++nXB; } 
          else            { std::cout<<"2"; ++nXA; }  
        }
        std::cout<<std::endl;
      }
      else{
        //sat -> print
        //std::cout<<status;
        assert(status==l_True|| nCi<=1);
        std::cout<<"PO "<<Abc_ObjName(pObj)<<" support partition: 0"<<std::endl;
      }
      

    }
    
    return 1;
}
