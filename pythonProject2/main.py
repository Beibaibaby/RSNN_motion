
structure=[]
def check (tmp,item):
    for i in tmp:
        if i == item: return False;

def get_on_table (structure):
 tmp =[]
 for item in structure:
   if first(item) == "On-Table" and (check (tmp,item)):
    tmp.append(item);


