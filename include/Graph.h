#ifndef GRAPH_H
#define GRAPH_H
#include <vector>
#include <algorithm>
using namespace std;

class VERTEX{
      public:
             // vertex id
             int id;
             // number of neighbor vertexes of this one
             int num_neighbor;
             // a erased vertex mark
             //bool is_deleted;
             // assign operator
             inline VERTEX& operator=(const VERTEX& to_assign)
             {
                    id = to_assign.id;
                    num_neighbor = to_assign.num_neighbor;
                    is_deleted_ = to_assign.is_deleted_;
                    return *this;
             }
             inline VERTEX():id(0),
                             num_neighbor(0),
                             is_deleted_(false){}
             inline VERTEX(const VERTEX& to_copy):id(to_copy.id),
                                                  num_neighbor(to_copy.num_neighbor),
                                                  is_deleted_(to_copy.is_deleted_){}
             bool is_deleted(){
                  return is_deleted_;
                  }
             void set_is_deleted(bool flag){
                  is_deleted_ = flag;
                  }
			 void clear(){
				 id = 0;
				 num_neighbor = 0;
				 is_deleted_ = false;
			 }
             private:
                     bool is_deleted_;
};

// to use sort algorithm in vector, define a utility method as the 3rd 
/*inline bool compare(const VERTEX& lhs, const VERTEX& rhs)
{
       return lhs.id < rhs.id;
}*/
// overload the "<" operator so as the VERTEX object can be sorted
inline bool operator<(const VERTEX& lhs, const VERTEX& rhs){
       return lhs.num_neighbor < rhs.num_neighbor;
       }
inline bool operator==(const VERTEX& lhs, const VERTEX& rhs){
                    return (lhs.id==rhs.id);
             }
             
class EDGE{
      public:
             // indicator of the first vertex
             int first;
             // indicator of the last vertex
             int last;
             // a vector holding the other vertex bypassed by this edge
			 vector<int> path;
             // a erased edge mark
             bool is_deleted;
			 EDGE():first(0),
				    last(0),
					path(),
					is_deleted(false){}
			 EDGE(const EDGE& edge):first(edge.first),
									last(edge.last),
									path(edge.path),
									is_deleted(edge.is_deleted){}
             bool has_vertex(VERTEX& v);
			 void clear(){
				 first = last = 0;
				 path.clear();
				 is_deleted = false;
			 }
};

inline bool operator==(const EDGE& lhs, const EDGE& rhs){
	return (lhs.first == rhs.first && lhs.last == rhs.last && lhs.path == rhs.path);
       }
inline bool EDGE::has_vertex(VERTEX& v){
     return (first == v.id || last == v.id);
     } 
inline EDGE contatenate(EDGE& edge1, EDGE& edge2, VERTEX& curr_vertex){
	EDGE new_edge;
	if (edge1.first == curr_vertex.id && edge1.last != curr_vertex.id)
	{
		swap(edge1.first, edge1.last);
	}
	if (edge2.last == curr_vertex.id && edge2.first != curr_vertex.id)
	{
		swap(edge2.first, edge2.last);
	}
	new_edge.first = edge1.first;
	new_edge.last = edge2.last;
	
	new_edge.path = edge1.path;
	new_edge.path.insert(new_edge.path.end(), edge2.path.begin(), edge2.path.end());
	new_edge.path.push_back(curr_vertex.id);
	/*sort(new_edge.path.begin(), new_edge.path.end());
	new_edge.path.erase(unique(new_edge.path.begin(), new_edge.path.end()), new_edge.path.end());*/
	new_edge.is_deleted = false;
	return new_edge;
}
#endif
             
