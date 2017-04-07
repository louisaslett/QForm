// Code to determine machine precision at a particular value.
// See https://en.wikipedia.org/wiki/Machine_epsilon#How_to_determine_machine_epsilon

typedef union {
  long long i64;
  double d64;
} dbl_64;

void machine_eps(double* value) {
  dbl_64 s;
  s.d64 = *value;
  s.i64++;
  *value = s.d64 - *value;
}
