namespace Signal {

  enum StopState {
    RUNNING,
    INTERRUPTED,
    STOPPING
  };

  enum Response {
    CONTINUE,
    DUMP,
    RESTART,
    STOP
  };

} // namespace Signal
