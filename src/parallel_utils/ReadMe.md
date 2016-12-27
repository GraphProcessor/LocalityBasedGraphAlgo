#Implementation Thoughts
##User Interfaces
- users are expected to use `map_reduce` or `map_reduce_fined_grained_redyce` or `map_fined_grained_reduce`, users are only 
required to use function template to get template functions

- `map_reduce` sample

```cpp
map_reduce(rand_access_user_collection, map_lambda, reduce_lambda)
```

- `map_fine_grained_reduce` sample, which requires reduced data to be collection of collection, e.g, `std::vector<std::vector<T>>`

```cpp
map_fine_grained_reduce(rand_acccess_user_collection, map_lambda, sucess_call_back_lambda, fail_call_back_lambda, flag_lambda)
```
