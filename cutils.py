

def sanify(*input):
    """ Makes sure that the name of the gwf target is not illegal. """
    input = ''.join([str(j) for j in input])

    output = []


    for i in input:
        
        ascii = ord(i)
        if (ascii >= 48 and ascii <= 57) or (ascii >= 65 and ascii <= 90) or (ascii >= 97 and ascii <= 122) or ascii == 95:
            output.append(i)
        else:
            output.append('_')


    return ''.join(output)


def dprint(*input):

    if not True:
        print(*input)
