def number_types(element_list):
    s, i = 0, 0
    for el in element_list:
        if el.type == 'S':
            s+=1
        elif el.type == 'I':
            i+=1

    return tuple([i, s])
