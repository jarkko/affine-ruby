class AffineTransformation
  def initialize(from, to)
    @from, @to = from, to

    if @from.size != @to.size || @from.size < 1
      raise "Both from and to should be of same size"
    end

    @dim = @to[0].size
    if @to.size < @dim
      raise "Too few points => under-determined system"
    end

    c = []
    (@dim + 1).times do
      c << Array.new(@dim, 0.0)
    end

    puts "c is #{c.inspect}"

    @dim.times do |j|
      (@dim + 1).times do |k|
        @from.size.times do |i|
          qt = @from[i] + [1]
          c[k][j] += qt[k] * @to[i][j]
        end
      end
    end

    q = []
    (@dim + 1).times do
      q << Array.new(@dim, 0.0) + [0]
    end

    @from.each do |qi|
      qt = qi + [1]
      (@dim + 1).times do |i|
        (@dim + 1).times do |j|
          q[i][j] += qt[i] * qt[j]
        end
      end
    end
    
    @m = []
    (@dim + 1).times do |i|
      @m << q[i] + c[i]
    end
    
    puts "@m is #{@m.inspect}"
    if !gauss_jordan(@m)
      raise "Singular matrix. Points are probably coplanar."
    end
  end
  
  def to_s
    res = ""
    @dim.times do |j|
      str = "x%d' = " % j
      @dim.times do |i|
        str += "x#{i} * #{@m[i][j+@dim+1]} + "
      end
      str += "%f" % @m[@dim][j + @dim + 1]
      res += str + "\n"
    end
    res
  end
  
  def transform(pt)
    res = Array.new(@dim, 0.0)
    @dim.times do |j|
      @dim.times do |i|
        res[j] += pt[i] * @m[i][j + @dim + 1]
      end
      res[j] += @m[@dim][j + @dim + 1]
    end
    res
  end

  private

  def gauss_jordan(m, eps = 1.0/(10**10))
    h, w = m.size, m[0].size
    h.times do |y|
      maxrow = y

      ((y+1)...h).each do |y2|
        if m[y2][y].abs > m[maxrow][y].abs
          maxrow = y2
        end
      end

      m[y], m[maxrow] = m[maxrow], m[y]
      if m[y][y].abs <= eps
        return false
      end

      ((y + 1)...h).each do |y2|
        c = m[y2][y] / m[y][y]
        (y...w).each do |x|
          m[y2][x] -= m[y][x] * c
        end
      end
    end

    y = h - 1
    while y >= 0
      c = m[y][y]
      (0...y).each do |y2|
        x = w - 1
        while x >= y
          m[y2][x] -= m[y][x] * m[y2][y] / c
          x -= 1
        end
      end
      m[y][y] /= c

      # normalize row y
      (h...w).each do |x|
        m[y][x] /= c
      end

      y -= 1
    end

    return true
  end
end

from = [[1,1], [1,2], [2,2], [2,1]]
to = [[4,4], [6,6], [8,4], [6,2]]

transformation = AffineTransformation.new(from, to)

puts "Transformation is #{transformation}"

err = 0.0

from.each_with_index do |i, idx|
  fp = i
  tp = to[idx]

  t = transformation.transform(fp)
  puts "#{fp.inspect} => #{t.inspect} ~= #{tp.inspect}"
  err += ((tp[0] - t[0])**2 + (tp[1] - t[1])**2)**0.5
end

puts "Fitting error = #{err}"