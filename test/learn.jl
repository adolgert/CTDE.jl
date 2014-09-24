# This explains how to make a base interface and derive from it.
# Not a class but an interface.

# Base module defines interface
module GSPN
export all_transitions
function all_transitions()
	return "base gspn"
end
end

# Derived module adds functions to methods of base.
module Simple
import GSPN.all_transitions
export A, all_transitions

type A
end

function all_transitions(x::A)
	"Simple A"
end

end


# Module using interface imports base, has access to derived.
module WithIt
using GSPN
function fire(a)
	println(all_transitions(a))
end
end


using Simple
using WithIt

a=Simple.A()
println(all_transitions(a))
WithIt.fire(a)