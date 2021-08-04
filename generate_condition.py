if __name__ == "__main__":
    conditions = open("base_cond.txt",mode="w")
    for i in range(1000):
        for j in range(1000):
            if(j == 0):
                conditions.write("400 ")
            elif(j == 999):
                conditions.write("100 ")
            else:
                conditions.write("200 ")
        conditions.write("\n")
    conditions.close()
