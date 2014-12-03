#!/usr/bin/env python

IGNORE_CHARACTERS = ["#", "-", "*", "^", "@"]

def removeIgnoreCharacters(string):
	for ig_char in IGNORE_CHARACTERS:
		string = string.replace(ig_char,"")
	return string

def averageList(l):
	return reduce(lambda x, y: x + y, l) / len(l)