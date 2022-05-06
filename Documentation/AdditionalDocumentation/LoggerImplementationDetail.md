Detailed explanation of the logger implementation by dducks for those interested.

A brief explanation how the logger class works. Beware, there is black magic going on here.

So, the previous version of the Logger used the syntax Logger<Log::LEVEL> logger; logger.log(Log::LEVEL, "Message",
args...); However, people didn't like the slightly large amount of log's in a single statement. Therefore, the
operator() has been chosen. The next problem I personally had with the logger is the fact that it relied on a
compile-time resolvable if-statement which could often (always) be optimised out. The problem however is that you would
need optimisation enabled to have no speed penalty for loglevels below visibility.

The problem however, is that we want to generate different code based on the first parameter; the loglevel. Because we
actually want to generate code and people don't like the preprocessor, we have to use the template system to do this.

One of the caveats of the template system is that it is only able to resolve templates based on parameters, not on
values. It also allows you to ommit template parameters in case where the template parameters can be deduced.

Based on these two rules, we now need not a value per loglevel, but a type. Therefore we use tagging, in which we use
the class LL with a template argument to create the different tags. Our loglevel, the enum class Log, is considered to
be a valid non-type template parameter - it must be either integral, enum
(which is integral) or a pointer - so we can use it to tag a class to create diffent types.

Now, our Logger instance has a template argument which is the loglevel associated with our Logger; anything of lesser
priority is ignored; It also honours the symbol MERCURY_LOGLEVEL, and compares it to the lesser of the template argument
and the preprocessor define.

The operator() function now should be based on the template parameter of the logger class (or MERCURY_LOGLEVEL), and the
loglevel of the current message. Because we can't depend on the value but just on the type, we have to do some magic.
Now, we want the function to resolve to something that produces output if the level is high enough, or nothing at all if
the level is of a too low priority for the current loglevel. For this we utilize std::enable_if<>
to select the right function.

Because in C++ a templated function which leads to ill-formed code upon instantiation is rejected and does NOT result in
a compile error, std::enable_if allows us to make either the version WITH output or WITHOUT output, based on template
parameters.

As you can see below, the class LL is completely empty; However, there are a few instances; one for every loglevel. Now,
remember, when using the logger, Log::DEFAULT resolves to the enum class value while DEFAULT resolves to the instance of
class LL with template argument Log::DEFAULT.

However, LL<Log::DEFAULT> differs in type from LL<Log::ERROR>, as the template argument is different. Because of
deduction, the compiler can figure out the actual values for the template arguments. In case of operator(), the first
argument is Log loglevel, which is deduced from the declaration;

template<Log LogLevel, typename... Args>
void operator(const LL<LogLevel> ll, std::string text, Args... args);

Now, the template parameter LogLevel gets filled with Log::ERROR in case we give ERROR as the first argument, or Log::
DEBUG with DEBUG as the first argument. Now, we want to resolve to two different functions, so as the return type we use
std::enable_if which would lead to ill-formed code in case the predicate of the std::enable_if is false. So, based on
the tag ll we now can select the different implementations.

It then buffers everything into a stringstream. so as long as operator<<
is defined properly for your object,. This then gets redirected towards the correct output channel.

Please note that operator() is an inline, templated function. In case the function body is empty, this code is very,
very likely to not emit any instructions at all. If you don't give arguments which have only non-const functions, the
function call can be considered invariant which means it can completely be taken out, in the case it's a seperate
function it just lowers the cost.
