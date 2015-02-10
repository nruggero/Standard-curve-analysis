#!/usr/bin/env python

IGNORE_CHARACTERS = ["#", "-", "*", "^", "@"]

def removeIgnoreCharacters(string):
	for ig_char in IGNORE_CHARACTERS:
		string = string.replace(ig_char,"")
	return string

def averageList(l):
	return reduce(lambda x, y: x + y, l) / len(l)

def findIndexClosestToNumber(number, myList):
	value_closest = min(myList, key=lambda x:abs(x-number))
	return myList.index(value_closest)

def median(mylist):
    sorts = sorted(mylist)
    length = len(sorts)
    if not length % 2:
        return (sorts[length / 2] + sorts[length / 2 - 1]) / 2.0
    return sorts[length / 2]