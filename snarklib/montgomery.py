"""Số học Montgomery: REDC, nhân Montgomery, chuyển vào/ra dạng Montgomery.

Với modulus lẻ p và R = 2^k > p (gcd(R, p) = 1), dạng Montgomery của a là
    ã = a·R mod p.
REDC(T) với 0 <= T < p·R trả về T·R^{-1} mod p chỉ bằng phép nhân, phép AND
với R - 1 (lấy mod R) và phép dịch phải k bit (chia cho R), không cần chia cho p:
    m = (T mod R)·p' mod R,   p' = -p^{-1} mod R
    t = (T + m·p) / R          (chia hết vì T + m·p ≡ 0 mod R)
    nếu t >= p: t = t - p       (bước trừ cuối, t < 2p)
Nhân Montgomery: mont_mul(ã, b̃) = REDC(ã·b̃) = (ab)·R mod p, tức là dạng
Montgomery của ab.

Hai lựa chọn R:
    Montgomery(p)            k = số bit của p (R nhỏ nhất là lũy thừa của 2 lớn hơn p)
    Montgomery(p, word=64)   k = 64·s, s = số từ máy 64 bit (cách arkworks làm:
                             BLS12-381 có p 381 bit -> s = 6, R = 2^384)
mul_cios cài phiên bản theo từng từ máy (CIOS, Koç–Acar–Kaliski 1996) và đếm
số phép nhân từ máy; kết quả phải trùng REDC(a·b).
"""


class Montgomery:
    def __init__(self, p, word=None):
        if p % 2 == 0 or p < 3:
            raise ValueError("modulus phải lẻ và >= 3")
        self.p = p
        bits = p.bit_length()
        if word is None:
            self.k = bits
            self.s = None
        else:
            self.s = -(-bits // word)          # số từ máy = ceil(bits / word)
            self.k = word * self.s
        self.word = word
        self.R = 1 << self.k
        self.mask = self.R - 1
        self.p_prime = (-pow(p, -1, self.R)) % self.R   # p' = -p^{-1} mod R
        self.R_mod_p = self.R % p
        self.R2 = (self.R * self.R) % p                 # dùng để chuyển vào dạng Montgomery
        self.so_lan_tru_cuoi = 0                        # đếm số lần bước trừ cuối xảy ra

    def __repr__(self):
        return "Montgomery(p = %d bit, R = 2^%d)" % (self.p.bit_length(), self.k)

    # --- REDC và nhân ---
    def redc(self, T):
        """T·R^{-1} mod p, yêu cầu 0 <= T < p·R."""
        assert 0 <= T < self.p * self.R
        m = ((T & self.mask) * self.p_prime) & self.mask
        t = (T + m * self.p) >> self.k
        if t >= self.p:
            t -= self.p
            self.so_lan_tru_cuoi += 1
        return t

    def mul(self, a_m, b_m):
        """Nhân hai phần tử dạng Montgomery: trả về (ab)R mod p."""
        return self.redc(a_m * b_m)

    def to_mont(self, a):
        """a -> a·R mod p, tính bằng REDC(a·R^2) (không chia cho p)."""
        return self.redc((a % self.p) * self.R2)

    def from_mont(self, a_m):
        """ã -> a = REDC(ã)."""
        return self.redc(a_m)

    def pow(self, a_m, e):
        """Lũy thừa trong dạng Montgomery (bình phương và nhân)."""
        r = self.R_mod_p                 # dạng Montgomery của 1
        for bit in bin(e)[2:]:
            r = self.mul(r, r)
            if bit == "1":
                r = self.mul(r, a_m)
        return r

    # --- CIOS theo từ máy ---
    def mul_cios(self, a_m, b_m):
        """Nhân Montgomery theo từng từ máy (CIOS). Trả về (kết quả, số phép nhân từ máy).

        Mỗi vòng i: t += a·b_i (s phép nhân), m = t_0·p'_0 mod W (1 phép nhân),
        t = (t + m·p) / W (s phép nhân) -> tổng 2s^2 + s phép nhân từ máy.
        """
        if self.word is None:
            raise ValueError("cần khởi tạo với word=...")
        w, s = self.word, self.s
        W = 1 << w
        wm = W - 1
        p0_prime = self.p_prime & wm        # p'_0 = -p^{-1} mod W
        a = [(a_m >> (w * j)) & wm for j in range(s)]
        pl = [(self.p >> (w * j)) & wm for j in range(s)]
        t = [0] * (s + 2)
        nhan = 0
        for i in range(s):
            bi = (b_m >> (w * i)) & wm
            C = 0
            for j in range(s):                       # t = t + a·b_i
                x = t[j] + a[j] * bi + C
                nhan += 1
                t[j], C = x & wm, x >> w
            x = t[s] + C
            t[s], t[s + 1] = x & wm, x >> w
            m = (t[0] * p0_prime) & wm               # m = t_0·p'_0 mod W
            nhan += 1
            x = t[0] + m * pl[0]
            nhan += 1
            C = x >> w
            for j in range(1, s):                    # t = (t + m·p) / W
                x = t[j] + m * pl[j] + C
                nhan += 1
                t[j - 1], C = x & wm, x >> w
            x = t[s] + C
            t[s - 1], C = x & wm, x >> w
            t[s] = t[s + 1] + C
        r = sum(t[j] << (w * j) for j in range(s + 1))
        if r >= self.p:
            r -= self.p
        return r, nhan
