////////////////////////////////////////////////////////////////////////////////
//! \file
//! \author Adam M Phillippy
//! \date 03/26/2003
//!
//! \brief Source for non-inline member functions of delta.hh
//!
//! \see delta.hh
////////////////////////////////////////////////////////////////////////////////

#include "delta.hh"
#include "fasta.hh"
#include "delcher.hh"
#include <cmath>
#include <sstream>
#include <algorithm>
using namespace std;



inline long int ScoreLocal
(long int scorej,
 long int leni, long int lenj,
 long int olap, float idyi, float maxolap);

  
struct LIS_t
//!< LIS score
{
  DeltaEdgelet_t * a;
  long int score;
  long int diff;
  long int from;
  bool used;
};


struct EdgeletQCmp_t
//!< Compares query lo coord
{
  bool operator( ) (const DeltaEdgelet_t * i, const DeltaEdgelet_t * j) const
  {
    //-- Sorting by score in the event of a tie ensures that when building
    //   LIS chains, the highest scoring ones get seen first, thus avoiding
    //   overlap problems

    if ( i -> loQ < j -> loQ )
      return true;
    else if ( i -> loQ > j -> loQ )
      return false;
    else if ( ScoreLocal (0, i->hiQ - i->loQ + 1, 0, 0, i->idy, 0) >
              ScoreLocal (0, j->hiQ - j->loQ + 1, 0, 0, j->idy, 0) )
      return true;
    else
      return false;
  }
};


struct EdgeletRCmp_t
//!< Compares reference lo coord
{
  bool operator( ) (const DeltaEdgelet_t * i, const DeltaEdgelet_t * j) const
  {
    //-- Sorting by score in the event of a tie ensures that when building
    //   LIS chains, the highest scoring ones get seen first, thus avoiding
    //   overlap problems

    if ( i -> loR < j -> loR )
      return true;
    else if ( i -> loR > j -> loR )
      return false;
    else if ( ScoreLocal (0, i->hiR - i->loR + 1, 0, 0, i->idy, 0) >
              ScoreLocal (0, j->hiR - j->loR + 1, 0, 0, j->idy, 0) )
      return true;
    else
      return false;
  }
};


struct NULLPred_t
//!< Return true if pointer is not NULL
{
  bool operator( ) (const void * i) const
  { return (i != NULL); }
};


//------------------------------------------------------------ DiffAligns ------
inline long int DiffAligns
(const DeltaEdgelet_t * i, const DeltaEdgelet_t * j)
{
  return
    ( ( j->loR < i->loR ) ? labs (j->hiR - i->loR) : labs (i->hiR - j->loR) ) +
    ( ( j->loQ < i->loQ ) ? labs (j->hiQ - i->loQ) : labs (i->hiQ - j->loQ) );
}


//------------------------------------------------------------ PickBest --------
long int PickBest
(const LIS_t * lis, const vector<long int> & allbest, float epsilon)
{
  long int size = allbest . size( );
  if ( epsilon < 0 && size != 0 )
    {
      long int eqc = 0;
      for ( ; eqc < size; ++ eqc )
        if ( lis[allbest[eqc]] . diff != lis[allbest . front( )] . diff )
          break;
      
      return (int)((double)eqc*rand() / (RAND_MAX + 1.0));
    }
  return size;
}


//------------------------------------------------------------------ RevC ------
inline unsigned long int RevC (const unsigned long int & coord,
                               const unsigned long int & len)
{
  return len - coord + 1;
}


//------------------------------------------------------------ ScoreLocal ------
inline long int ScoreLocal
(long int scorej, long int leni, long int lenj,
long int olap, float idyi, float maxolap)
{
  if ( olap > 0  &&
       ((float)olap / (float)leni * 100.0 > maxolap  ||
	(float)olap / (float)lenj * 100.0 > maxolap) )
    return -1;
  else
    return (scorej + (long int)((leni - olap) * pow (idyi, 2)));
}


//----------------------------------------------------------- ScoreGlobal ------
inline long int ScoreGlobal
(long int scorej, long int leni, long int olap, float idyi)
{
  return (scorej + (long int)((leni - olap) * pow (idyi, 2)));
}


//------------------------------------------------------------------ Swap ------
inline void Swap (unsigned long int & a, unsigned long int & b)
{
  unsigned long int t = a; a = b; b = t;
}


//------------------------------------------------------------ UpdateBest ------
bool UpdateBest
(LIS_t * lis, long int size, vector<long int> & allbest, float epsilon)
{
  if ( size == 0 )
    return false;

  long int best, i;

  //-- Find the best
  for ( best = 0; best < size; ++ best )
    if ( ! lis[best] . used )
      break;
  for ( i = best + 1; i < size; ++ i )
    if ( ! lis[i] . used
         &&
         ( lis[i] . score > lis[best] . score
           ||
           (lis[i] . score == lis[best] . score &&
            lis[i] . diff  <  lis[best] . diff) ) )
      best = i;

  //-- Nonequivalent
  if ( ! allbest . empty( )
       &&
       (best == size
        ||
        (epsilon < 0 &&
         lis[best].score < lis[allbest.front( )].score)
        ||
        (epsilon >= 0 &&
         (float)(lis[allbest . front( )].score - lis[best].score) /
         (float)(lis[allbest . front( )].score) * 100.0 > epsilon)) )
    return false;
  
  //-- Equivalent
  allbest . push_back (best);

  for ( i = best; i >= 0  &&  i < size; i = lis[i] . from )
    lis[i] . used = true;

  return true;
}


//===================================================== DeltaReader_t ==========
//----------------------------------------------------- open -------------------
void DeltaReader_t::open
     (const string & delta_path)
{
  delta_path_m = delta_path;

  //-- Open the delta file
  delta_stream_m.open (delta_path_m.c_str ( ));
  checkStream ( );

  //-- Read the file header
  delta_stream_m >> reference_path_m;
  delta_stream_m >> query_path_m;
  delta_stream_m >> data_type_m;
  if ( (data_type_m != NUCMER_STRING  &&  data_type_m != PROMER_STRING) )
    delta_stream_m.setstate (ios::badbit);
  checkStream ( );
  is_open_m = true;

  //-- Advance to first record header
  while ( delta_stream_m.peek ( ) != '>' )
    if ( delta_stream_m.get ( ) == EOF )
      break;
}


//----------------------------------------------------- readNextAlignment ------
void DeltaReader_t::readNextAlignment
     (DeltaAlignment_t & align, const bool read_deltas)
{
  long int delta;
  float total;

  //-- Make way for the new alignment
  align.clear ( );

  //-- Read the alignment header
  delta_stream_m >> align.sR;
  delta_stream_m >> align.eR;
  delta_stream_m >> align.sQ;
  delta_stream_m >> align.eQ;
  delta_stream_m >> align.idyc;
  delta_stream_m >> align.simc;
  delta_stream_m >> align.stpc;
  if ( align.sR <= 0  ||  align.eR <= 0  ||
       align.sQ <= 0  ||  align.eQ <= 0  ||
       align.idyc < 0  ||  align.simc < 0  ||  align.stpc < 0 )
    delta_stream_m.setstate (ios::badbit);
  checkStream ( );

  total = labs(align.eR - align.sR) + 1.0;
  if ( data_type_m == PROMER_STRING )
    total /= 3.0;

  //-- Get all the deltas
  do
    {
      delta_stream_m >> delta;
      checkStream ( );

      if ( delta < 0 )
	total ++;
      if ( read_deltas )
	align.deltas.push_back (delta);
    } while ( delta != 0 );

  //-- Flush the remaining whitespace
  while ( delta_stream_m.get ( ) != '\n' );

  //-- Calculate the identity, similarity and stopity
  align.idy = (total - (float)align.idyc) / total * 100.0;
  align.sim = (total - (float)align.simc) / total * 100.0;
  align.stp = (float)align.stpc / (total * 2.0) * 100.0;
}


//----------------------------------------------------- readNextRecord ---------
bool DeltaReader_t::readNextRecord (const bool read_deltas)
{
  //-- If EOF or any other abnormality
  if ( delta_stream_m.peek ( ) != '>' )
    return false;

  //-- Make way for the new record
  record_m.clear ( );
  is_record_m = true;

  //-- Read the record header
  delta_stream_m.get ( );
  delta_stream_m >> record_m.idR;
  delta_stream_m >> record_m.idQ;
  delta_stream_m >> record_m.lenR;
  delta_stream_m >> record_m.lenQ;
  if ( record_m.lenR <= 0  ||  record_m.lenQ <= 0 )
    delta_stream_m.setstate (ios::badbit);
  checkStream ( );

  //-- Flush the remaining whitespace
  while ( delta_stream_m.get ( ) != '\n' );

  //-- For each alignment...
  DeltaAlignment_t align;
  while ( delta_stream_m.peek ( ) != '>'  &&
	  delta_stream_m.peek ( ) != EOF )
    {
      readNextAlignment (align, read_deltas);
      record_m.aligns.push_back (align);
    }

  return true;
}


//===================================================== DeltaEdge_t ============
//------------------------------------------------------build ------------------
void DeltaEdge_t::build (const DeltaRecord_t & rec)
{
  stringstream ss;
  vector<long int>::const_iterator di;
  DeltaEdgelet_t * p;

  vector<DeltaAlignment_t>::const_iterator i;
  for ( i = rec . aligns . begin( ); i != rec . aligns . end( ); ++ i )
    {
      //-- Set the edgelet
      p = new DeltaEdgelet_t( );

      p -> edge = this;

      p -> idy = i -> idy / 100.0;
      p -> sim = i -> sim / 100.0;
      p -> stp = i -> stp / 100.0;

      p -> idyc = i -> idyc;
      p -> simc = i -> simc;
      p -> stpc = i -> stpc;

      p -> loR = i -> sR;
      p -> hiR = i -> eR;
      p -> loQ = i -> sQ;
      p -> hiQ = i -> eQ;

      p -> dirR = p -> hiR < p -> loR ? REVERSE_DIR : FORWARD_DIR;
      p -> dirQ = p -> hiQ < p -> loQ ? REVERSE_DIR : FORWARD_DIR;

      //-- Get the delta information
      for ( di = i -> deltas . begin( ); di != i -> deltas . end( ); ++ di )
        {
          ss << *di << '\n';
          p -> delta . append (ss . str( ));
          ss . str ("");
        }

      //-- Force loR < hiR && loQ < hiQ
      if ( p -> dirR == REVERSE_DIR )
        Swap (p -> loR, p -> hiR);
      if ( p -> dirQ == REVERSE_DIR )
        Swap (p -> loQ, p -> hiQ);

      edgelets . push_back (p);
    }
}


//===================================================== DeltaGraph_t ===========
//----------------------------------------------------- build ------------------
void DeltaGraph_t::build (const string & deltapath, bool getdeltas)
{
  DeltaReader_t dr;
  DeltaEdge_t * dep;
  pair<map<string, DeltaNode_t>::iterator, bool> insret;


  //-- Open the delta file and read in the alignment information
  dr . open (deltapath);

  refpath = dr . getReferencePath( );
  qrypath = dr . getQueryPath( );

  if ( dr . getDataType( ) == NUCMER_STRING )
    datatype = NUCMER_DATA;
  else if ( dr . getDataType( ) == PROMER_STRING )
    datatype = PROMER_DATA;
  else
    datatype = NULL_DATA;

  //-- Read in the next graph edge, i.e. a new delta record
  while ( dr . readNext (getdeltas) )
    {
      dep = new DeltaEdge_t( );


      //-- Find the reference node in the graph, add a new one if necessary
      insret = refnodes . insert
        (map<string, DeltaNode_t>::value_type
         (dr.getRecord( ).idR, DeltaNode_t( )));
      dep -> refnode = &((insret . first) -> second);

      //-- If a new reference node
      if ( insret . second )
        {
          dep -> refnode -> id  = &((insret . first) -> first);
          dep -> refnode -> len = dr . getRecord( ) . lenR;
        }


      //-- Find the query node in the graph, add a new one if necessary
      insret = qrynodes . insert
        (map<string, DeltaNode_t>::value_type
         (dr.getRecord( ).idQ, DeltaNode_t( )));
      dep -> qrynode = &((insret . first) -> second);

      //-- If a new query node
      if ( insret . second )
        {
          dep -> qrynode -> id  = &((insret . first) -> first);
          dep -> qrynode -> len = dr . getRecord( ) . lenQ;
        }


      //-- Build the edge
      dep -> build (dr . getRecord( ));
      dep -> refnode -> edges . push_back (dep);
      dep -> qrynode -> edges . push_back (dep);
    }
  dr . close ( );
}


//------------------------------------------------------------------- clean ----
void DeltaGraph_t::clean( )
{
  map<string, DeltaNode_t>::iterator i;
  map<string, DeltaNode_t>::iterator ii;
  vector<DeltaEdge_t *>::iterator j;
  vector<DeltaEdgelet_t *>::iterator k;

  //-- For all ref nodes
  for ( i = refnodes . begin( ); i != refnodes . end( ); )
    {
      //-- For all edges
      for ( j  = i -> second . edges . begin( );
            j != i -> second . edges . end( ); ++ j )
        {
          //-- For all edgelets
          for ( k  = (*j) -> edgelets . begin( );
                k != (*j) -> edgelets . end( ); ++ k )
            {
              if ( ! (*k) -> isGOOD )
                {
                  //-- Delete the bad edgelet and mark for erasure
                  delete (*k);
                  *k = NULL;
                }
            }

          //-- Erase the marked edgelets
          k = partition ((*j) -> edgelets . begin( ),
                         (*j) -> edgelets . end( ),
                         NULLPred_t( ));
          (*j) -> edgelets . erase (k, (*j) -> edgelets . end( ));

          //-- Mark the edge if empty
          if ( (*j) -> edgelets . empty( ) )
            *j = NULL;
        }

      //-- Erase the marked edges
      j = partition (i -> second . edges . begin( ),
                     i -> second . edges . end( ),
                     NULLPred_t( ));
      i -> second . edges . erase (j, i -> second . edges . end( ));

      //-- Erase the node if empty
      ii = i ++;
      if ( ii -> second . edges . empty( ) )
        refnodes . erase (ii);
    }

  //-- For all qry nodes
  for ( i = qrynodes . begin( ); i != qrynodes . end( ); )
    {
      for ( j  = i -> second . edges . begin( );
            j != i -> second . edges . end( ); ++ j )
        {
          //-- Delete the edge if empty and mark for erasure
          if ( (*j) -> edgelets . empty( ) )
            {
              delete (*j);
              *j = NULL;
            }
        }

      //-- Erase the marked edges
      j = partition (i -> second . edges . begin( ),
                     i -> second . edges . end( ),
                     NULLPred_t( ));
      i -> second . edges . erase (j, i -> second . edges . end( ));

      //-- Erase the node if empty
      ii = i ++;
      if ( ii -> second . edges . empty( ) )
        qrynodes . erase (ii);
    }
 
}


//------------------------------------------------------------ getNodeCount ----
long int DeltaGraph_t::getNodeCount( )
{
  long int sum = refnodes . size( ) + qrynodes . size( );
  return sum;
}


//------------------------------------------------------------ getEdgeCount ----
long int DeltaGraph_t::getEdgeCount( )
{
  long int sum = 0;

  map<string, DeltaNode_t>::iterator i;
  for ( i = refnodes . begin( ); i != refnodes . end( ); ++ i )
    sum += i -> second . edges . size( );

  return sum;
}


//--------------------------------------------------------- getEdgeletCount ----
long int DeltaGraph_t::getEdgeletCount( )
{
  long int sum = 0;

  map<string, DeltaNode_t>::iterator i;
  vector<DeltaEdge_t *>::iterator j;
  for ( i = refnodes . begin( ); i != refnodes . end( ); ++ i )
    for ( j  = i -> second . edges . begin( );
          j != i -> second . edges . end( ); ++ j )
      sum += (*j) -> edgelets . size( );

  return sum;
}


//-------------------------------------------------------------- flagGLIS ------
void DeltaGraph_t::flagGLIS (float epsilon)
{
  LIS_t * lis = NULL;
  long int lis_size = 0;
  long int i, j, n;
  long int olap, olapQ, olapR, len, lenQ, lenR, score, diff;

  vector<DeltaEdgelet_t *> edgelets;

  map<string, DeltaNode_t>::iterator mi;
  vector<DeltaEdge_t *>::iterator ei;
  vector<DeltaEdgelet_t *>::iterator eli;


  //-- For each reference sequence
  for ( mi = refnodes . begin( ); mi != refnodes.end( ); ++ mi )
    {
      //-- For each query aligning to this reference
      for ( ei  = (mi -> second) . edges . begin( );
            ei != (mi -> second) . edges . end( ); ++ ei )
        {
          //-- Clear any previous chains
          (*ei) -> chains . clear( );

          //-- Collect all the good edgelets
          edgelets . clear( );
          for ( eli  = (*ei) -> edgelets . begin( );
                eli != (*ei) -> edgelets . end( ); ++ eli )
            {
              //-- Clear any previous chain info
              (*eli) -> next = NULL;

              if ( (*eli) -> isGOOD )
                {
                  edgelets . push_back (*eli);
                  
                  //-- Fix the coordinates to make global LIS work
                  if ( (*eli) -> dirR == (*eli) -> dirQ )
                    {
                      (*eli) -> dirQ = FORWARD_DIR;
                    }
                  else
                    {
                      if ( (*eli) -> dirQ == REVERSE_DIR )
                        Swap ((*eli) -> loQ, (*eli) -> hiQ);
                      (*eli) -> loQ = RevC ((*eli) -> loQ, (*ei)->qrynode->len);
                      (*eli) -> hiQ = RevC ((*eli) -> hiQ, (*ei)->qrynode->len);
                      (*eli) -> dirQ = REVERSE_DIR;
                    }
                }
            }

          //-- Resize and initialize
          n = edgelets . size( );
          if ( n > lis_size )
            {
              lis = (LIS_t *) Safe_realloc (lis, sizeof (LIS_t) * n);
              lis_size = n;
            }
          for ( i = 0; i < n; ++ i )
            lis[i] . used = false;

          //-- Sort by lo query coord
          sort (edgelets . begin( ), edgelets . end( ), EdgeletQCmp_t( ));

          //-- Continue until all equivalent repeats are extracted
          vector<long int> allbest;
          do
            {
              //-- Dynamic
              for ( i = 0; i < n; ++ i )
                {
                  if ( lis[i] . used ) continue;

                  lis[i] . a = edgelets[i];
                  
                  lenR = lis[i] . a -> hiR - lis[i] . a -> loR + 1;
                  lenQ = lis[i] . a -> hiQ - lis[i] . a -> loQ + 1;
                  len = lenR > lenQ ? lenQ : lenR;
                  lis[i] . score = ScoreGlobal (0, len, 0, lis[i] . a -> idy);

                  lis[i] . from = -1;
                  lis[i] . diff = 0;
                  
                  for ( j = 0; j < i; ++ j )
                    {
                      if ( lis[j] . used ) continue;

                      if ( lis[i] . a -> dirQ != lis[j] . a -> dirQ )
                        continue;
                  
                      lenR = lis[i] . a -> hiR - lis[i] . a -> loR + 1;
                      lenQ = lis[i] . a -> hiQ - lis[i] . a -> loQ + 1;
                      len = lenR > lenQ ? lenQ : lenR;
                  
                      olapR = lis[j] . a -> hiR - lis[i] . a -> loR + 1;
                      olapQ = lis[j] . a -> hiQ - lis[i] . a -> loQ + 1;
                      olap = olapR > olapQ ? olapR : olapQ;
                      if ( olap < 0 )
                        olap = 0;

                      diff = lis[j] . diff + DiffAligns (lis[i].a, lis[j].a);

                      score = ScoreGlobal
                        (lis[j] . score, len, olap, lis[i] . a -> idy);

                      if ( score > lis[i] . score
                           ||
                           (score == lis[i] . score && diff < lis[i] . diff) )
                        {
                          lis[i] . from = j;
                          lis[i] . score = score;
                          lis[i] . diff = diff;
                        }
                    }
                }
            } while ( UpdateBest (lis, n, allbest, epsilon) );

          long int beg = PickBest (lis, allbest, epsilon);
          long int end = allbest . size( );
          if ( beg == end ) beg = 0;
          else end = beg + 1;

          //-- Flag the edgelets
          DeltaEdgelet_t * last;
          for ( ; beg < end; ++ beg )
            {
              last = NULL;
              for ( i = allbest[beg]; i >= 0  &&  i < n; i = lis[i] . from )
                {
                  lis[i] . a -> isGLIS = true;
                  lis[i] . a -> next = last;
                  last = lis[i] . a;
                }
              (*ei) -> chains . push_back (last);
            }

          //-- Repair the coordinates
          for ( eli = edgelets . begin( ); eli != edgelets . end( ); ++ eli )
            {
              if ( ! (*eli) -> isGLIS )
                (*eli) -> isGOOD = false;

              if ( (*eli) -> dirQ == FORWARD_DIR )
                {
                  (*eli) -> dirQ = (*eli) -> dirR;
                }
              else
                {
                  if ( (*eli) -> dirR == FORWARD_DIR )
                    Swap ((*eli) -> loQ, (*eli) -> hiQ);
                  (*eli) -> loQ = RevC ((*eli) -> loQ, (*ei)->qrynode->len);
                  (*eli) -> hiQ = RevC ((*eli) -> hiQ, (*ei)->qrynode->len);
                  (*eli) -> dirQ =
                    (*eli) -> dirR == FORWARD_DIR ? REVERSE_DIR : FORWARD_DIR;
                }
            }
        }
    }

  free (lis);
}


//------------------------------------------------------------- flagScore ------
void DeltaGraph_t::flagScore (long int minlen, float minidy)
{
  map<string, DeltaNode_t>::const_iterator mi;
  vector<DeltaEdge_t *>::const_iterator ei;
  vector<DeltaEdgelet_t *>::iterator eli;

  for ( mi = refnodes . begin( ); mi != refnodes . end( ); ++ mi )
    for ( ei  = (mi -> second) . edges . begin( );
          ei != (mi -> second) . edges . end( ); ++ ei )
      for ( eli  = (*ei) -> edgelets . begin( );
            eli != (*ei) -> edgelets . end( ); ++ eli )
        if ( (*eli) -> isGOOD )
          {
            //-- Flag low identities
            if ( (*eli)->idy * 100.0 < minidy )
              (*eli) -> isGOOD = false;

            //-- Flag small lengths
            if ( (*eli)->hiR - (*eli)->loR + 1 < (unsigned long)minlen ||
                 (*eli)->hiQ - (*eli)->loQ + 1 < (unsigned long)minlen )
              (*eli) -> isGOOD = false;
          }
}


//-------------------------------------------------------------- flagQLIS ------
void DeltaGraph_t::flagQLIS (float epsilon, float maxolap)
{
  LIS_t * lis = NULL;
  long int lis_size = 0;
  long int i, j, n;
  long int olap, leni, lenj, score, diff;

  vector<DeltaEdgelet_t *> edgelets;

  map<string, DeltaNode_t>::iterator mi;
  vector<DeltaEdge_t *>::iterator ei;
  vector<DeltaEdgelet_t *>::iterator eli;


  //-- For each query sequence
  for ( mi = qrynodes . begin( ); mi != qrynodes . end( ); ++ mi )
    {
      //-- Clean any previous chains
      (mi -> second) . chains . clear( );

      //-- For each reference aligning to this query
      edgelets . clear( );
      for ( ei  = (mi -> second) . edges . begin( );
            ei != (mi -> second) . edges . end( ); ++ ei )
        {
          //-- Collect all the good edgelets
          for ( eli  = (*ei) -> edgelets . begin( );
                eli != (*ei) -> edgelets . end( ); ++ eli )
            {
              //-- Clean any previous chain info
              (*eli) -> next = NULL;
              
              if ( (*eli) -> isGOOD )
                edgelets . push_back (*eli);
            }
        }

      //-- Resize and initialize
      n = edgelets . size( );
      if ( n > lis_size )
        {
          lis = (LIS_t *) Safe_realloc (lis, sizeof (LIS_t) * n);
          lis_size = n;
        }
      for ( i = 0; i < n; ++ i )
        lis[i] . used = false;

      //-- Sort by lo query coord
      sort (edgelets . begin( ), edgelets . end( ), EdgeletQCmp_t( ));

      //-- Continue until all equivalent repeats are extracted
      vector<long int> allbest;
      do
        {
          //-- Dynamic
          for ( i = 0; i < n; ++ i )
            {
              if ( lis[i] . used ) continue;

              lis[i] . a = edgelets[i];

              leni = lis[i] . a -> hiQ - lis[i] . a -> loQ + 1;
              lis[i] . score =
                ScoreLocal (0, leni, 0, 0, lis[i].a->idy, 0);

              lis[i] . from = -1;
              lis[i] . diff = 0;

              for ( j = 0; j < i; ++ j )
                {
                  if ( lis[j] . used ) continue;

                  leni = lis[i] . a -> hiQ - lis[i] . a -> loQ + 1;
                  lenj = lis[j] . a -> hiQ - lis[j] . a -> loQ + 1;
                  olap = lis[j] . a -> hiQ - lis[i] . a -> loQ + 1;
                  if ( olap < 0 )
                    olap = 0;

                  diff = lis[j] . diff + DiffAligns (lis[i].a, lis[j].a);

                  score = ScoreLocal
                    (lis[j] . score, leni, lenj, olap, lis[i].a->idy, maxolap);

                  if ( score > lis[i] . score
                       ||
                       (score == lis[i] . score && diff < lis[i] . diff) )
                    {
                      lis[i] . from = j;
                      lis[i] . score = score;
                      lis[i] . diff = diff;
                    }
                }
            }
        } while ( UpdateBest (lis, n, allbest, epsilon) );

      long int beg = PickBest (lis, allbest, epsilon);
      long int end = allbest . size( );
      if ( beg == end ) beg = 0;
      else end = beg + 1;

      //-- Flag the edgelets
      DeltaEdgelet_t * last;
      for ( ; beg < end; ++ beg )
        {
          last = NULL;
          for ( i = allbest[beg]; i >= 0  &&  i < n; i = lis[i] . from )
            {
              lis[i] . a -> isQLIS = true;
              lis[i] . a -> next = last;
              last = lis[i] . a;
            }
          (mi -> second) . chains . push_back (last);
        }

      for ( eli = edgelets . begin( ); eli != edgelets . end( ); ++ eli )
        if ( ! (*eli) -> isQLIS )
          (*eli) -> isGOOD = false;
    }

  free (lis);
}


//-------------------------------------------------------------- flagRLIS ------
void DeltaGraph_t::flagRLIS (float epsilon, float maxolap)
{
  LIS_t * lis = NULL;
  long int lis_size = 0;
  long int i, j, n;
  long int olap, leni, lenj, score, diff;

  vector<DeltaEdgelet_t *> edgelets;

  map<string, DeltaNode_t>::iterator mi;
  vector<DeltaEdge_t *>::iterator ei;
  vector<DeltaEdgelet_t *>::iterator eli;


  //-- For each reference sequence
  for ( mi = refnodes . begin( ); mi != refnodes . end( ); ++ mi )
    {
      //-- Clean any previous chains
      (mi -> second) . chains . clear( );

      //-- For each query aligning to this reference
      edgelets . clear( );
      for ( ei  = (mi -> second) . edges . begin( );
            ei != (mi -> second) . edges . end( ); ++ ei )
        {
          //-- Collect all the good edgelets
          for ( eli  = (*ei) -> edgelets . begin( );
                eli != (*ei) -> edgelets . end( ); ++ eli )
            {
              //-- Clean any previous chain info
              (*eli) -> next = NULL;

              if ( (*eli) -> isGOOD )
                edgelets . push_back (*eli);
            }
        }

      //-- Resize
      n = edgelets . size( );
      if ( n > lis_size )
        {
          lis = (LIS_t *) Safe_realloc (lis, sizeof (LIS_t) * n);
          lis_size = n;
        }
      for ( i = 0; i < n; ++ i )
        lis[i] . used = false;

      //-- Sort by lo reference coord
      sort (edgelets . begin( ), edgelets . end( ), EdgeletRCmp_t( ));

      //-- Continue until all equivalent repeats are extracted
      vector<long int> allbest;
      do
        {
          //-- Dynamic
          for ( i = 0; i < n; ++ i )
            {
              if ( lis[i] . used ) continue;

              lis[i] . a = edgelets[i];

              leni = lis[i] . a -> hiR - lis[i] . a -> loR + 1;
              lis[i] . score =
                ScoreLocal (0, leni, 0, 0, lis[i].a->idy, 0);

              lis[i] . from = -1;
              lis[i] . diff = 0;

              for ( j = 0; j < i; ++ j )
                {
                  if ( lis[j] . used ) continue;

                  leni = lis[i] . a -> hiR - lis[i] . a -> loR + 1;
                  lenj = lis[j] . a -> hiR - lis[j] . a -> loR + 1;
                  olap = lis[j] . a -> hiR - lis[i] . a -> loR + 1;
                  if ( olap < 0 )
                    olap = 0;

                  diff = lis[j] . diff + DiffAligns (lis[i].a, lis[j].a);

                  score = ScoreLocal
                    (lis[j] . score, leni, lenj, olap, lis[i].a->idy, maxolap);

                  if ( score > lis[i] . score
                       ||
                       (score == lis[i] . score && diff < lis[i] . diff) )
                    {
                      lis[i] . from = j;
                      lis[i] . score = score;
                      lis[i] . diff = diff;
                    }
                }
            }
        } while ( UpdateBest (lis, n, allbest, epsilon) );

      long int beg = PickBest (lis, allbest, epsilon);
      long int end = allbest . size( );
      if ( beg == end ) beg = 0;
      else end = beg + 1;
      
      //-- Flag the edgelets
      DeltaEdgelet_t * last;
      for ( ; beg < end; ++ beg )
        {
          last = NULL;
          for ( i = allbest[beg]; i >= 0  &&  i < n; i = lis[i] . from )
           {
             lis[i] . a -> isRLIS = true;
             lis[i] . a -> next = last;
             last = lis[i] . a;
           }
          (mi -> second) . chains . push_back (last);
        }

      for ( eli = edgelets . begin( ); eli != edgelets . end( ); ++ eli )
        if ( ! (*eli) -> isRLIS )
          (*eli) -> isGOOD = false;
    }

  free (lis);
}


//-------------------------------------------------------------- flagUNIQ ------
void DeltaGraph_t::flagUNIQ (float minuniq)
{
  unsigned long int i, uniq, len;

  vector<DeltaEdgelet_t *> edgelets;

  map<string, DeltaNode_t>::const_iterator mi;
  vector<DeltaEdge_t *>::const_iterator ei;
  vector<DeltaEdgelet_t *>::iterator eli;


  //-- For each reference sequence
  unsigned long int ref_size = 0;
  unsigned long int ref_len = 0;
  unsigned char * ref_cov = NULL;
  for ( mi = refnodes . begin( ); mi != refnodes . end( ); ++ mi )
    {
      //-- Reset the reference coverage array
      ref_len = (mi -> second) . len;
      if ( ref_len > ref_size )
        {
          ref_cov = (unsigned char *) Safe_realloc (ref_cov, ref_len + 1);
          ref_size = ref_len;
        }
      for ( i = 1; i <= ref_len; ++ i )
        ref_cov[i] = 0;

      //-- Collect all the good edgelets
      edgelets . clear( );
      for ( ei  = (mi -> second) . edges . begin( );
            ei != (mi -> second) . edges . end( ); ++ ei )
        for ( eli  = (*ei) -> edgelets . begin( );
              eli != (*ei) -> edgelets . end( ); ++ eli )
          if ( (*eli) -> isGOOD )
            {
              edgelets . push_back (*eli);

              //-- Add to the reference coverage
              for ( i = (*eli) -> loR; i <= (*eli) -> hiR; i ++ )
                if ( ref_cov[i] < UCHAR_MAX )
                  ref_cov[i] ++;
            }

      //-- Calculate the uniqueness of each edgelet
      for ( eli = edgelets . begin( ); eli != edgelets . end( ); ++ eli )
        {
          uniq = 0;
          len = (*eli) -> hiR - (*eli) -> loR + 1;
          for ( i = (*eli) -> loR; i <= (*eli) -> hiR; i ++ )
            if ( ref_cov[i] == 1 )
              uniq ++;

          //-- Flag low reference uniqueness
          if ( (float)uniq / (float)len * 100.0 < minuniq )
            (*eli) -> isGOOD = false;
        }
    }
  free (ref_cov);


  //-- For each query sequence
  unsigned long int qry_size = 0;
  unsigned long int qry_len = 0;
  unsigned char * qry_cov = NULL;
  for ( mi = qrynodes . begin( ); mi != qrynodes . end( ); ++ mi )
    {
      //-- Reset the query coverage array
      qry_len = (mi -> second) . len;
      if ( qry_len > qry_size )
        {
          qry_cov = (unsigned char *) Safe_realloc (qry_cov, qry_len + 1);
          qry_size = qry_len;
        }
      for ( i = 1; i <= qry_len; ++ i )
        qry_cov[i] = 0;

      //-- Collect all the good edgelets
      edgelets . clear( );
      for ( ei  = (mi -> second) . edges . begin( );
            ei != (mi -> second) . edges . end( ); ++ ei )
        for ( eli  = (*ei) -> edgelets . begin( );
              eli != (*ei) -> edgelets . end( ); ++ eli )
          if ( (*eli) -> isGOOD )
            {
              edgelets . push_back (*eli);

              //-- Add to the query coverage
              for ( i = (*eli) -> loQ; i <= (*eli) -> hiQ; i ++ )
                if ( qry_cov[i] < UCHAR_MAX )
                  qry_cov[i] ++;
            }

      //-- Calculate the uniqueness of each edgelet
      for ( eli = edgelets . begin( ); eli != edgelets . end( ); ++ eli )
        {
          uniq = 0;
          len = (*eli) -> hiQ - (*eli) -> loQ + 1;
          for ( i = (*eli) -> loQ; i <= (*eli) -> hiQ; i ++ )
            if ( qry_cov[i] == 1 )
              uniq ++;
          
          //-- Flag low query uniqueness
          if ( (float)uniq / (float)len * 100.0 < minuniq )
            (*eli) -> isGOOD = false;
        }
    }
  free (qry_cov);
}


//----------------------------------------------------- loadSequences ----------
void DeltaGraph_t::loadSequences ( )
{
  const int init_size = 10000;
  const int max_line = 1024;

  map<string, DeltaNode_t>::iterator mi;

  FILE * qryfile, * reffile;
  char * R = NULL;
  char * Q = NULL;
  char id [max_line];
  long int initsize;
  long int len;

  //-- Read in the reference sequences
  reffile = File_Open (refpath . c_str( ), "r");
  initsize = init_size;
  R = (char *) Safe_malloc (initsize);
  while ( Fasta_Read_String (reffile, R, initsize, id, FALSE) )
    if ( (mi = refnodes . find (id)) != refnodes . end( ) )
      {
        len = strlen (R + 1);
        mi -> second . seq = (char *) Safe_malloc (len + 2);
        mi -> second . seq[0] = '\0';
        strcpy (mi -> second . seq + 1, R + 1);
        if ( (unsigned long) len != mi -> second . len )
          {
            cerr << "ERROR: Reference input does not match delta file\n";
            exit (EXIT_FAILURE);
          }
      }
  fclose (reffile);
  free (R);

  //-- Read in the query sequences
  qryfile = File_Open (qrypath . c_str( ), "r");
  initsize = init_size;
  Q = (char *) Safe_malloc (initsize);
  while ( Fasta_Read_String (qryfile, Q, initsize, id, FALSE) )
    if ( (mi = qrynodes . find (id)) != qrynodes . end( ) )
      {
        len = strlen (Q + 1);
        mi -> second . seq = (char *) Safe_malloc (len + 2);
        mi -> second . seq[0] = '\0';
        strcpy (mi -> second . seq + 1, Q + 1);
        if ( (unsigned long) len != mi -> second . len )
          {
            cerr << "ERROR: Query input does not match delta file\n";
            exit (EXIT_FAILURE);
          }
      }
  fclose (qryfile);
  free (Q);


  //-- Check that we found all the sequences
  for ( mi = refnodes.begin( ); mi != refnodes.end( ); ++ mi )
    if ( mi -> second . seq == NULL )
      {
        cerr << "ERROR: '" << mi -> first << "' not found in reference file\n";
        exit (EXIT_FAILURE);
      }

  for ( mi = qrynodes.begin( ); mi != qrynodes.end( ); ++ mi )
    if ( mi -> second . seq == NULL )
      {
        cerr << "ERROR: '" << mi -> first << "' not found in query file\n";
        exit (EXIT_FAILURE);
      }
}


//----------------------------------------------------- outputDelta ------------
ostream & DeltaGraph_t::outputDelta (ostream & out)
{
  bool header;
  unsigned long int s1, e1, s2, e2;
 
  map<string, DeltaNode_t>::const_iterator mi;
  vector<DeltaEdge_t *>::const_iterator ei;
  vector<DeltaEdgelet_t *>::const_iterator eli;
 
  //-- Print the file header
  cout
    << refpath << ' ' << qrypath << '\n'
    << (datatype == PROMER_DATA ? PROMER_STRING : NUCMER_STRING) << '\n';
 
  for ( mi = qrynodes . begin( ); mi != qrynodes . end( ); ++ mi )
    {
      for ( ei  = (mi -> second) . edges . begin( );
            ei != (mi -> second) . edges . end( ); ++ ei )
        {
          header = false;
 
          for ( eli  = (*ei) -> edgelets . begin( );
                eli != (*ei) -> edgelets . end( ); ++ eli )
            {
              if ( ! (*eli) -> isGOOD )
                continue;
 
              //-- Print the sequence header
              if ( ! header )
                {
                  cout
                    << '>'
                    << *((*ei) -> refnode -> id) << ' '
                    << *((*ei) -> qrynode -> id) << ' '
                    << (*ei) -> refnode -> len << ' '
                    << (*ei) -> qrynode -> len << '\n';
                  header = true;
                }
              //-- Print the alignment
              s1 = (*eli) -> loR;
              e1 = (*eli) -> hiR;
              s2 = (*eli) -> loQ;
              e2 = (*eli) -> hiQ;
              if ( (*eli) -> dirR == REVERSE_DIR )
                Swap (s1, e1);
              if ( (*eli) -> dirQ == REVERSE_DIR )
                Swap (s2, e2);

              cout
                << s1 << ' ' << e1 << ' ' << s2 << ' ' << e2 << ' '
                << (*eli) -> idyc << ' '
                << (*eli) -> simc << ' '
                << (*eli) -> stpc << '\n'
                << (*eli) -> delta;
            }
        }
    }
  return out;
}
