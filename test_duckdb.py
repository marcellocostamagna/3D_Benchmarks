import duckdb

con = duckdb.connect(":memory:")
con.execute("CREATE TABLE X(i INT)")
data = [(x,) for x in range(1, 10_000)]
con.execute("INSERT INTO X VALUES (?)", data)
