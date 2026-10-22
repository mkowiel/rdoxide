"""Compatibility types shared by the Boost.Python and the nanobind flavours of
the RDKit Python API.

RDKit has two wrapper sets whose tests contradict each other in a few places:
the Boost.Python wrappers return tuples where the nanobind ones return lists,
and they raise different exception types for the same misuse. rdoxide
satisfies both with these types: a tuple that also compares equal to the
equivalent list, and exceptions deriving from both exception classes.
"""


class TupleList(tuple):
  """A tuple that compares equal to a list (or tuple) with equal items."""

  __slots__ = ()

  def __eq__(self, other):
    if isinstance(other, (tuple, list)):
      return len(self) == len(other) and all(a == b for a, b in zip(self, other))
    return NotImplemented

  def __ne__(self, other):
    r = self.__eq__(other)
    return r if r is NotImplemented else not r

  __hash__ = tuple.__hash__


class OverflowTypeError(OverflowError, TypeError):
  """Boost.Python raises OverflowError, nanobind TypeError."""


class ValueTypeError(ValueError, TypeError):
  """Boost.Python raises ValueError, nanobind TypeError."""


class RuntimeTypeError(RuntimeError, TypeError):
  """Boost.Python raises TypeError, nanobind RuntimeError."""
