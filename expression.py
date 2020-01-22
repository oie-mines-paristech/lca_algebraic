
from typing import Dict, Union, List

class Expression() :

    # To be overridden
    def __operator__(self, operand, other):
        if not isinstance(other, Expression) :
            if isinstance(other, float) or isinstance(other, int):
                other = Litteral(other)
            else:
                raise Exception("Operand not supported : %s (%s)" % (other, type(other)))

        return OperatorExpression(operand, [self, other])

    def __add__(self, other):
        return self.__operator__('+', other)

    def __sub__(self, other):
        return self.__operator__('-', other)

    def __mul__(self, other):
        return self.__operator__('*', other)

    def __truediv__(self, other):
        return self.__operator__('/', other)

    def __pow__(self, power):
        return self.__operator__('^', power)

    def __neg__(self):
        return Negative(self)



class Litteral(Expression) :

    def __init__(self, value):
        self.value = value

    def __repr__(self):
        return str(self.value)

class Negative(Expression) :
    def __init__(self, operand):
        self.operand = operand
    def __repr__(self):
        return "-%s" % repr(self.operand)

class OperatorExpression(Expression) :

    def __operator__(self, operand, other):
        # In case of commutative operand, append arguments
        if operand == self.operand and operand in ['+', '*'] :
            self.args.append(other)
            return self
        else :
            return super(OperatorExpression, self).__operator__(operand, other)

    def __init__(self, operand:str, args:List[Expression]):
        self.args = args
        self.operand = operand

    def __repr__(self):
        args = ["(%s)" % repr(arg) if isinstance(arg, OperatorExpression) else repr(arg) for arg in self.args]
        return (" %s " % self.operand).join(args)
