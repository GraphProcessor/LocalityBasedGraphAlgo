#Thoughts on Parallel Util
##Things to Clarify
- where to put the data, how to extract data

##Implementation Issues
        http://drdobbs.com/cpp/232500059
        
        There are, indeed, performance issues with std:function that must be taken into account whenever using it. The main strength of std::function, namely, its type-erasure mechanism, does not come for free, and we might (but not necessarily must) pay a price for that.
        
        std::function is a template class that wraps callable types. However, it is not parametrized on the callable type itself but only on its return and argument types. The callable type is known only at construction time and, therefore, std::function cannot have a pre-declared member of this type to hold a copy of the object given to its constructor.
        
        Roughly speaking (actually, things are more complicated than that) std::function can hold only a pointer to the object passed to its constructor, and this raises a lifetime issue. If the pointer points to an object whose lifetime is smaller than that of the std::function object, then the inner pointer will become dangling. To prevent this problem std::function might make a copy of the object on the heap through a call to operator new (or a custom allocator). The dynamic memory allocation is what people refer the most as a performance penalty implied by std::function.
        
        I have recently written an article with more details and that explains how (and where) one can avoid paying the price of a memory allocation.

