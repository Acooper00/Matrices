#! /usr/bin/python3

# Class for matricies

import math

class MatrixError(Exception):
    pass

class Matrix():

    def __init__(self,list_):

        self.matrix = []# This reassignement of the list is to prevent quasi-connectivity
        for row in list_:
            newrow = []
            for value in row:
                newrow.append(value)
            self.matrix.append(newrow)

        if self.matrix not in [[], [[]]]:
            self.size = [len(self.matrix[0]),len(self.matrix)]
        else:
            self.matrix = []
            self.size = [0,0]

        for row in range(0,len(self.matrix)):
            if type(self.matrix[row]) != type([]):
                raise MatrixError
            for value in self.matrix[row]:
                if type(value) not in [type(0),type(complex(0)),type(0.0)]:
                    raise MatrixError
            if len(self.matrix[row]) != self.size[0]:
                raise MatrixError

        return

    def __repr__(self):
        
        string = ""
        for row in self.matrix:
            string += "("
            for value in row:
                string += str(value).center(5," ")
            string += ")\n"
        string = string[:-1]

        return string

    def __add__(self,obj):# Overwriting default methods as pointers in order to give the ability to use operators -- see "https://realpython.com/operator-function-overloading/"
        return self.Add(obj)

    def __sub__(self,obj):
        return self.Subtract(obj)

    def __mul__(self,obj):
        return self.Multiply(obj)

    def __rmul__(self,obj):
        if type(obj) not in [type(0),type(complex(0)),type(0.0)]:
            raise MatrixError
        else:
            return self.Multiply(obj)

    def __pow__(self,obj):
        
        try:
            if int(obj) == -1:
                return self.Inverse()
            elif int(obj) >= 1:
                newmatrix = self.Identity()
                for i in range(0,int(obj)):
                    newmatrix = newmatrix.Multiply(self)
                return newmatrix
            else:
                raise MatrixError
        except ValueError:
            if str(obj).upper() == "T":
                return self.Transpose()
            else:
                return MatrixError

    def Value(self,row,collumn):

        return self.matrix[row][collumn]

    def Values(self):

        return list(self.matrix)

    def Dimensions(self):

        return str(self.size[1])+" x "+str(self.size[0])

    def Alter(self,row,collumn,value):

        if type(value) not in [type(0),type(complex(0)),type(0.0)]:
            raise MatrixError

        self.matrix[row][collumn] = value

        return

    def DeleteRow(self,row):

        del self.matrix[row]
        
        self.size[1] -= 1

        return

    def DeleteCollumn(self,collumn):

        table = list(self.matrix)

        for row in table:
            del row[collumn]

        self.matrix = list(table)

        self.size[0] -= 1

        return

    def Delete(self,type_,number):

        if type_.lower() == "row":
            self.DeleteRow(number)
        elif type_.lower() == "collumn":
            self.DeleteCollumn(number)
        else:
            raise MatrixError

        return

    def AddRow(self,row="",position=""):

        if row == "":
            if self.size[0] == 0:
                raise MatrixError
            row = []
            for i in range(self.size[0]):
                row.append(0)

        if position == "":
            position = self.size[1]

        if type(position) not in [type(0),type(complex(0)),type(0.0)]:
            raise MatrixError

        if position < 0 or position > self.size[1]:
            raise MatrixError

        if self.size[0] != 0:
            if len(row) != self.size[0]:
                raise MatrixError

        for value in row:
            if type(value) not in [type(0),type(complex(0)),type(0.0)]:
                raise MatrixError

        self.matrix.insert(position,list(row))

        if self.size[0] == 0:
            self.size[0] = len(row)
        self.size[1] += 1

        return

    def AddCollumn(self,collumn="",position=""):

        if collumn == "":
            if self.size[1] == 0:
                raise MatrixError
            collumn = []
            for i in range(self.size[1]):
                collumn.append(0)

        if position == "":
            position = self.size[0]

        if type(position) not in [type(0),type(complex(0)),type(0.0)]:
            raise MatrixError

        if position < 0 or position > self.size[0]:
            raise MatrixError

        if self.size[1] != 0:
            if len(collumn) != self.size[1]:
                raise MatrixError

        for value in collumn:
            if type(value) not in [type(0),type(complex(0)),type(0.0)]:
                raise MatrixError

        if self.size[1] == 0:
            for i in range(0,len(collumn)):
                self.matrix.append([])
        for row,value in zip(range(0,len(collumn)),collumn):
            self.matrix[row].insert(position,value)

        if self.size[1] == 0:
            self.size[1] = len(collumn)
        self.size[0] += 1

        return
    
    def Transpose(self):

        self.matrix = list(self.Transposition().Values())

        return

    def IsSquare(self):

        if self.size[0] == self.size[1]:
            return True
        elif self.size[0] != self.size[1]:
            return False

    def Zero(self):

        row = []
        for y in range(0,self.size[0]):
            row.append(0)

        zero = []
        for x in range(0,self.size[1]):
            zero.append(list(row))

        return Matrix(zero)

    def Identity(self):

        if self.IsSquare():
            identity = list(self.Zero().Values())
            for x in range(0,self.size[1]):
                identity[x][x] = 1

            return Matrix(identity)

        else:
            raise MatrixError

    def Transposition(self):
        
        newmatrix = Matrix([])

        for row in self.matrix:
            newmatrix.AddCollumn(row)

        return newmatrix

    def Cofactor(self,row,collumn):

        minor = self.Minor(row,collumn)
        sign = (-1)**(row+collumn)

        return sign*minor

    def Determinant(self):
        
        if self.IsSquare():

            if self.size[0] == 1:
                return self.Value(0,0)

            else:
                total = 0
                for collumn in range(0,self.size[1]):
                    total += self.Value(0,collumn)*self.Cofactor(0,collumn)

                return total
        else:
            return MatrixError

    def Minor(self,row,collumn):


        minor = Matrix(list(self.Values()))
        minor.DeleteRow(row)
        minor.DeleteCollumn(collumn)

        determinant = minor.Determinant()

        return determinant

    def Inverse(self):

        if self.IsSquare() != True:
            return MatrixError

        if self.Determinant() == 0:
            raise MatrixError

        newmatrix = Matrix(list(self.matrix))

        det = self.Determinant()
        for row in range(0,newmatrix.size[1]):
            for collumn in range(0,newmatrix.size[0]):
                newmatrix.matrix[row][collumn] = self.Cofactor(row,collumn)/det

        newmatrix.Transpose()

        return newmatrix

    def Add(self,matrix):

        if type(matrix) != type(self):
            raise MatrixError

        if self.size != matrix.size:
            raise MatrixError

        newmatrix = Matrix(self.matrix)

        for row in range(0,self.size[1]):
            for collumn in range(0,self.size[0]):
                newmatrix.matrix[row][collumn] += matrix.matrix[row][collumn]

        return newmatrix

    def Subtract(self,matrix):

        if type(matrix) != type(self):
            raise MatrixError

        if self.size != matrix.size:
            raise MatrixError

        newmatrix = Matrix(self.matrix)

        for row in range(0,self.size[1]):
            for collumn in range(0,self.size[0]):
                newmatrix.matrix[row][collumn] -= matrix.matrix[row][collumn]

        return newmatrix

    def Multiply(self,factor):

        if type(factor) == type(self):

            if self.size[0] != factor.size[1]:
                raise MatrixError

            newmatrix = Matrix([[]])
            for row in range(0,self.size[1]):
                newmatrix.AddRow([0 for i in range(newmatrix.size[0])])
            for collumn in range(0,factor.size[0]):
                newmatrix.AddCollumn([0 for i in range(newmatrix.size[1])])

            for row in range(0,self.size[1]):
                for collumn in range(0,factor.size[0]):
                    for i in range(0,self.size[0]):
                        newmatrix.matrix[row][collumn] += self.matrix[row][i]*factor.matrix[i][collumn]

            return newmatrix

        elif type(factor) in [type(0),type(complex(0)),type(0.0)]:

            newmatrix = self.Zero()

            for row in range(0,self.size[1]):
                for collumn in range(0,self.size[0]):
                    newmatrix.matrix[row][collumn] = self.matrix[row][collumn]*factor

            return newmatrix

        else:
            raise MatrixError

# A Set Of Common Matrices For Usage

UnitX = Matrix([[1],[0]])
UnitY = Matrix([[0],[1]])

PauliX = Matrix([[0,1],[1,0]])
PauliY = Matrix([[0,complex(0,-1)],[complex(0,1),0]])
PauliZ = Matrix([[1,0],[0,-1]])

Hadamard = (1/math.sqrt(2))*Matrix([[1,1],[1,-1]])

def Rotation2D(angle, unit="rad"):

    if unit == "rad":
        pass
    elif unit == "deg":
        angle = (math.pi/180)*angle
    else:
        raise ValueError

    rotation = Matrix([[math.cos(angle), (-1)*math.sin(angle)],[math.sin(angle), math.cos(angle)]])

    return rotation
