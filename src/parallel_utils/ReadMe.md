#Parallel Utils
##Static Schedule, Two Stage Execution
1. fixed number of seeds are transformed into intermediate communities, which are going to be merged into a global one. 
tasks are split statically. intermediate representation should be able to be iterated in community-as-an-element manner. 

2. merging a community into a community list is regarded as a bulk-synchronization execution, where the number of community 
list will grow with computation. detailed execution context is a thread pool with synchronization and callback support.
