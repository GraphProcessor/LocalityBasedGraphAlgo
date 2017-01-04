#Design Thoughts
##User Interfaces
- users are expected to use `map_reduce` or `map_reduce_fined_grained_redyce` or `map_fined_grained_reduce`, users are only 
required to use function template to get template functions

- `map_reduce` sample

```cpp
map_reduce(rand_access_user_collection, map_lambda, reduce_lambda)
```

- `map_fine_grained_reduce` sample, which requires reduced data to be collection of collection, e.g, `std::vector<std::vector<T>>`

```cpp
map_fine_grained_reduce(rand_acccess_user_collection, map_lambda, sucess_call_back_lambda, fail_call_back_lambda, predicate_lambda)
```

##Implementation
###Map-Reduce
- make use of directed data graph in `map_reduce` to assign tasks, seeking for data locality and load balance.
- the scheduling of reduce tasks could be based on certain comparison criteria, e.g, sort by collection size if the reduce task is in Container Type.

###Map-Reduce with Parallel STL
- make use of map-reduce data-flow and parallel stl for `set_intersection`

###Map-FIne-Grained-Reduce
- since the type to be reduce here is not a simple value type, instead it is a collection type, the size of which may increase with the 
progress of reduce phase, e.g, `{{1,2,3}}, ...` at first with small size, but `{{...},...}` at last with huge number of communities and 
each community has huge number of members

- make use of thread pool for fine-grained-reduce, which makes tasks a little bit balanced in last few comparisons, where computations are
 only set_intersection operation, due to the nature of locality-based overlapping community detection feature

