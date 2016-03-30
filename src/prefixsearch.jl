"""
This is a binary tree where the leaves are values
and the nodes are sums of those values. It is meant to make it
easier to find the leaf such that the sum of it and all previous
leaves is greater than a given value.
"""
type PrefixSearchTree{T<:Real}
	array::Array{T,1}
	depth::Int
	offset::Int
end

"""
Constructor for an empty tree of type ElementType and length count.
"""
function PrefixSearchTree(ElementType::DataType, count::Int)
	PrefixSearchTree(zeros(ElementType, count))
end

"""
Constructor of a prefix search tree from an iterable list of real numbers.
"""
function PrefixSearchTree(summables)
	a=Array(summables)
	assert(length(a)>0)
	depth=Int(ceil(log(2, length(a))))+1
	offset=2^(depth-1)
	b=zeros(eltype(summables), 2^depth - 1)
	b[offset:(length(a)+offset-1)]=a
	#print("depth $depth offset $offset length $(length(b))\n")
	pst=PrefixSearchTree{eltype(summables)}(b, depth, offset)
	CalculatePrefix(pst)
	pst
end

"""
Find the minimum index such that the prefix is greater than the given
value.

Precondition: The value must be strictly less than the total for
the tree.
"""
function Choose(pst::PrefixSearchTree, value)
	if value>=pst.array[1]
		error("Value $value not less than total $(pst.array[1])")
	end

	requested_value=value
    index=1
    for level = 0:(pst.depth-2)
        left_child=2*index
        if pst.array[left_child]>value
            index=left_child
        else
            index=left_child+1
            value-=pst.array[left_child]
        end
    end
    index-pst.offset+1, pst.array[index]-value
end

Total(pst::PrefixSearchTree)=pst.array[1]

"""
If there are multiple values to enter, then present them
at once as pairs of tuples, (index, value).
"""
function Update!(pst::PrefixSearchTree, pairs)
	modify=Set{Int}()
	for (pos, value) in pairs
		index=pos+pst.offset-1
		pst.array[index]=value
		push!(modify, div(index, 2))
	end

	# everything at depth-1 is correct, and changes are in modify.
    for depth = range(pst.depth-2, -1, pst.depth-1)
        parents=Set{Int}()
        for node_idx in modify
            pst.array[node_idx]=pst.array[2*node_idx]+pst.array[2*node_idx+1]
            push!(parents, div(node_idx, 2))
        end
        modify=parents
    end
end

function Update!{T}(pst::PrefixSearchTree{T}, index::Int, value::T)
    Update!(pst, [(index, value)])
end

# Private
function CalculatePrefix(pst::PrefixSearchTree)
	for d = range(pst.depth-2, -1, pst.depth-1)
		for node_idx = (2^d):(2^(d+1)-1)
			pst.array[node_idx]=pst.array[2*node_idx]+pst.array[2*node_idx+1]
		end
	end
end


function TestSingleValue()
	t=PrefixSearchTree([3])
	assert(Total(t)==3)
	assert(Choose(t, 2)[1]==1)
	Update!(t, [(1,4)])
	assert(Total(t)==4)
	assert(Choose(t, 3.3)[1]==1)
end

function TestTwoValues()
	t=PrefixSearchTree([3, 1])
	assert(Total(t)==4)
	assert(Choose(t, 2)[1]==1)
	assert(Choose(t, 3)[1]==2)
	assert(Choose(t, 3.7)[1]==2)
	Update!(t, [(1, 4)])
	v=[(2,1), (3.3,1), (4.1,2)]
	for (guess, result) in v
		assert(Choose(t, guess)[1]==result)
	end
	Update!(t, [(2,2)])
	v=[(2, 1), (3.3, 1), (5.1, 2)]
	for (guess, result) in v
		assert(Choose(t, guess)[1]==result)
	end
end


function TestThreeValues()
    t=PrefixSearchTree([3, 1, 2])
    assert(Total(t) == 6)
    assert(Choose(t, 2.1)[1] == 1)
    assert(Choose(t, 3)[1] == 2)
    assert(Choose(t, 3.1)[1] == 2)
    assert(Choose(t, 5)[1] == 3)
    Update!(t, ((2,2),(3,3)))
    v=[(2,1),(3,2),(4.8,2),(5.1,3),(7.9,3)]
    for (guess, result) in v
        assert(Choose(t, guess)[1] == result)
    end
end

function TestThreeFloats()
    t=PrefixSearchTree([3.5, 1.5, 2.5])
    assert(abs(Total(t)-7.5)<1e-9)
    assert(Choose(t, 2.1)[1] == 1)
    assert(Choose(t, 3.5)[1] == 2)
    assert(Choose(t, 3.7)[1] == 2)
    assert(Choose(t, 5.5)[1] == 3)
    Update!(t, ((2,2.5),(3,3.5)))
    v=[(2,1),(3.5,2),(4.8,2),(6.1,3),(7.9,3)]
    for (guess, result) in v
        assert(Choose(t, guess)[1] == result)
    end
end

function TestFourValues()
    vals=[3,1,2,4]
    t=PrefixSearchTree(vals)
    assert(Total(t) == 10)
    v=[(2,1),(3,2),(4.1,3),(9.9,4)]
    for (guess, result) in v
        assert(Choose(t, guess)[1] == result)
    end
end


function TestFiveValues()
    vals=[3,0,2,4,1]
    t=PrefixSearchTree(vals)
    assert(Total(t) == 10)
    v=[(2,1),(3,3),(4.1,3),(9.9,5)]
    for (guess, result) in v
        assert(Choose(t, guess)[1] == result)
    end
    Update!(t, [(2,7),(3,3)])
    assert(Total(t) == 18)
    v=[(17.7,5)]
    for (guess, result) in v
        assert(Choose(t, guess)[1] == result)
    end
end

function TestPrefixSearchTree()
	TestSingleValue()
	TestTwoValues()
	TestThreeValues()
	TestThreeFloats()
	TestFourValues()
	TestFiveValues()
end
