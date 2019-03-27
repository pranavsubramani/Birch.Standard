/**
 * Event handler that replays a trace of events in immediate mode.
 *
 * - replay: The trqce to replay. This may be either an immediate trace or a
 *   delayed trace. For the latter, any eliminated variables will be
 *   simulated as needed.
 */
class ImmediateReplayHandler(replay:List<Event>) < ReplayHandler(replay) {
  function handle(evt:FactorEvent) -> Real {
    return evt.observe();
  }
  
  function handle(evt:RandomEvent) -> Real {
    replayEvt:Event?;
    if !replay.empty() {
      replayEvt <- replay.front();
      replay.popFront();
    }    
    if evt.hasValue() {
      ///@todo Check that values match
      return evt.observe();
    } else if replayEvt? {
      evt.value(replayEvt!);
    } else {
      evt.value();
    }
    return 0.0;
  }
}
